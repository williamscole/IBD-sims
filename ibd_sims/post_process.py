import argparse
import importlib
import os
import yaml
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Callable, Optional
import shutil
import submitit

from utils import apply_overrides

# ── Sub-configs ───────────────────────────────────────────────────────────────
class PostProcessConfig:
    
    def __init__(self, name):
        self.dir_name = name

class AnalysisConfig:

    # Top-level resource defaults
    _DEFAULTS = {
        "local": True,
        "workers": 1,
        "mem_gb": 8,
        "time_min": 60,
    }

    def __init__(self, args: dict, path: str):

        # Set resource defaults (overwritten if present in args)
        for key, val in self._DEFAULTS.items():
            setattr(self, key, val)

        # Gets the analyses to do
        if args.get("post_process") is None or args["post_process"] == "none":
            analysis_names = []
        else:
            analysis_names = [a.strip() for a in args["post_process"].split(",")]

        for top_key, top_val in args.items():
            if isinstance(top_val, dict) and top_key in analysis_names:
                pp_obj = PostProcessConfig(top_key)
                
                for key, val in top_val.items():
                    setattr(pp_obj, key, val)

                setattr(self, top_key, pp_obj)

            else:
                setattr(self, top_key, top_val)

        self.analyses = analysis_names
        self.path = path

    @classmethod
    def from_yaml(cls, yaml_path, overrides=None):

        with open(str(yaml_path), 'r') as file:
            args = yaml.safe_load(file)

        if overrides:
            args = apply_overrides(args, overrides)

        args["n_iter"] = args["iter"]
        del args["iter"]

        return cls(args, os.path.dirname(str(yaml_path)))

    def dump_config(self, out_path):
        """Dump the current analysis config (with overrides) to a YAML file."""
        d = {}
        for key in vars(self):
            val = getattr(self, key)
            if isinstance(val, PostProcessConfig):
                d[key] = {k: v for k, v in vars(val).items()}
            elif key not in ("analyses", "path"):
                d[key] = val
        # Restore iter key name for consistency
        if "n_iter" in d:
            d["iter"] = d.pop("n_iter")
        with open(out_path, "w") as f:
            yaml.dump(d, f, default_flow_style=False)

# ── PostProcessor ABC ─────────────────────────────────────────────────────────

class PostProcessor(ABC):
    """Base class for analysis steps that run on completed simulation output.

    Subclasses must set `sub_config_key` to the name of their sub-config
    attribute on AnalysisConfig (e.g., "ibdne"). This is used to resolve
    resource settings: sub-config values override top-level defaults.
    """

    sub_config_key: str  # Subclasses must define this

    def __init__(self, config: AnalysisConfig, path: str, n_iter: int = None, iter_n: int = None):
        self.config = config
        self.path = path

        # Exactly one must not be None
        assert (n_iter is None) + (iter_n is None) == 1
    
        self.n_iter = n_iter
        self.iter_n = iter_n

        # If true, just a single iteration needs to be run; as given by iter_n
        self.single_iter = n_iter is None

    def _get_sub_config(self):
        """Return the sub-config dataclass for this processor."""
        return getattr(self.config, self.sub_config_key)

    def _get_resource(self, name: str):
        """Resolve a resource value: sub-config overrides top-level.

        Checks self.config.<sub_config_key>.<name> first; if None,
        falls back to self.config.<name>.
        """
        sub = self._get_sub_config()
        value = getattr(sub, name, None)
        if value is not None:
            return value
        return getattr(self.config, name)

    # Fields to ignore when comparing configs for reuse.
    # Subclasses should override this with their resource-only fields.
    resource_fields: list = []

    def make_dir(self, path, dir_name):
        made_dir = False
        out_dir = ""
        cfg = self._get_sub_config()

        for i in range(1, 1000):
            target_dir = os.path.join(path, dir_name, f"{i:03d}")

            if os.path.exists(target_dir):
                yaml_file = os.path.join(target_dir, "args.yaml")

                if not os.path.exists(yaml_file):
                    continue

                tmp_config = AnalysisConfig.from_yaml(yaml_file)
                tmp_sub = getattr(tmp_config, self.sub_config_key, None)

                if tmp_sub is None:
                    continue

                # Compare all non-resource attributes
                is_match = True
                for attr in vars(cfg):
                    if attr in self.resource_fields:
                        continue
                    if getattr(cfg, attr, None) != getattr(tmp_sub, attr, None):
                        is_match = False
                        break

                if is_match:
                    out_dir = target_dir
                    made_dir = True
                    break
            else:
                os.makedirs(target_dir, exist_ok=True)
                out_dir = target_dir
                made_dir = True
                break

        assert made_dir, f"Could not create directory: all slots 001-999 are full in {path}/{dir_name}"

        self.out_dir = out_dir

        return out_dir

    def _execute_helper(self):
        cfg = self._get_sub_config()
        out_dir = self.make_dir(self.path, cfg.dir_name)
        self.config.dump_config(f"{out_dir}/args.yaml")
        return out_dir

    def _submit_jobs(self):
        """Submit Slurm jobs and return them without waiting."""
        executor = submitit.SlurmExecutor(folder=f"{self.path}/slurm")
        executor.update_parameters(
            mem=f"{self._get_resource('mem_gb')}GB",
            time=self._get_resource("time_min"),
            cpus_per_task=self._get_resource("workers"),
            use_srun=False,
        )
        iters = list(range(1, self.n_iter + 1))
        return executor.map_array(self._single_iter, iters)

    def _execute_loop(self, wait=True):
        local = self._get_resource("local")
        iters = list(range(1, self.n_iter + 1))

        if local:
            workers = self._get_resource("workers")
            if workers > 1:
                from concurrent.futures import ProcessPoolExecutor
                with ProcessPoolExecutor(max_workers=workers) as pool:
                    list(pool.map(self._single_iter, iters))
            else:
                for iter_n in iters:
                    self._single_iter(iter_n)
        else:
            jobs = self._submit_jobs()
            if wait:
                for job in jobs:
                    job.result()
            else:
                for job in jobs:
                    print(f"Submitted job {job.job_id}")

    def _single_iter(self, iter_n):
        """Run the analysis for a single iteration. Override in subclass."""
        ...

    @abstractmethod
    def execute(self, wait=True):
        """Run the analysis."""
        ...

# ── Orchestrator ──────────────────────────────────────────────────────────────

def postprocess(config_or_args, n_iter=None, iter_n=None, path=None, wait=True):
    if isinstance(config_or_args, AnalysisConfig):
        config = config_or_args
    else:
        config = AnalysisConfig(config_or_args, path)

    processors = []
    for key in config.analyses:
        module_name = getattr(config, key).path.replace(".py", "")
        class_name = getattr(config, key).object
        mod = importlib.import_module(module_name)
        cls = getattr(mod, class_name)
        processor = cls(config, path, n_iter=n_iter, iter_n=iter_n)
        processors.append(processor)

    local = any(p._get_resource("local") for p in processors)

    if local:
        for p in processors:
            p.execute(wait=wait)
    else:
        # Submit all analyses concurrently
        all_jobs = []
        for p in processors:
            p._execute_helper()  # create output dirs first
            all_jobs.extend(p._submit_jobs())

        if wait:
            for job in all_jobs:
                job.result()
        else:
            for job in all_jobs:
                print(f"Submitted job {job.job_id}")


# ── CLI helpers ───────────────────────────────────────────────────────────────

def _coerce(value: str):
    """Cast a CLI string value to the appropriate Python type."""
    if value.lower() == "true":
        return True
    if value.lower() == "false":
        return False
    if value.lower() in ("null", "none"):
        return None
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    return value


def main():
    parser = argparse.ArgumentParser(
        description="Run post-processing analyses on a completed simulation."
    )
    parser.add_argument(
        "input",
        help="Path to an analysis YAML config, or a run directory containing args.yaml",
    )
    parser.add_argument(
        "--n-iter", type=int, default=None,
        help="Number of iterations (default: read from args.yaml in run directory)",
    )
    parser.add_argument(
        "--no-wait", action="store_true", default=False,
        help="Submit Slurm jobs and exit immediately without waiting for them to finish.",
    )
    parser.add_argument(
        "--set", nargs="*", metavar="KEY=VALUE", default=None, action="append",
        help=(
            "Override config values. Use dot notation for nested keys. "
            "Examples: --set ibdne.nboots=100  --set local=false  --set ibdne.workers=4"
        ),
    )

    args = parser.parse_args()

    # Determine if input is a YAML file or a run directory
    input_path = Path(args.input)

    if input_path.is_dir() and (input_path / "args.yaml").exists():
        yaml_file = input_path / "args.yaml"

    elif input_path.is_file():
        yaml_file = input_path

    else:
        parser.error(f"Input not found: {args.input}")

    config = AnalysisConfig.from_yaml(yaml_file, overrides=args.set)

    postprocess(config, n_iter=config.n_iter, path=config.path, wait=not args.no_wait)


if __name__ == "__main__":
    main()