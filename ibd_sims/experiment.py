import yaml
import sys
import os
from pathlib import Path


# ── YAML I/O ──────────────────────────────────────────────────────────────────

def load_exp_yaml(yaml_file):
    exp_dict = yaml.safe_load(open(yaml_file))
    return exp_dict


def write_out_yaml(shared_args, exp_args, resources, yaml_path):
    args_dict = shared_args.copy()

    # Merge resources (max rule)
    for resource, resource_val in resources.items():
        args_dict[resource] = max(resource_val, args_dict.get(resource, 0))

    # Merge exp_args, stripping any nested 'resources' blocks
    for i, j in exp_args.items():
        args_dict[i] = {k: v for k, v in j.items() if k != "resources"} if isinstance(j, dict) else j

    with open(yaml_path, "w") as outf:
        yaml.dump(args_dict, outf, default_flow_style=False)

    return args_dict


# ── Label ─────────────────────────────────────────────────────────────────────

def get_label(demo_name, mating_name, end_chr=None, all_end_chrs=None):
    label = f"{demo_name}__{mating_name}"
    if end_chr is not None and all_end_chrs is not None and len(all_end_chrs) > 1:
        label += f"__chr{end_chr}"
    return label


# ── Resource extraction ───────────────────────────────────────────────────────

def extract_resources(*args):
    out_resources = {}
    for arg in args:
        if isinstance(arg, dict) and "resources" in arg:
            for resource, resource_val in arg["resources"].items():
                out_resources[resource] = max(out_resources.get(resource, 0), resource_val)
    return out_resources


# ── Setup ─────────────────────────────────────────────────────────────────────

def setup_dirs(exp_dict):
    base_dir = Path(exp_dict["experiment"])
    yaml_dir = base_dir / "yaml_files"
    yaml_dir.mkdir(parents=True, exist_ok=True)
    exp_dict["base_dir"] = str(base_dir)
    return yaml_dir


# ── YAML generation ───────────────────────────────────────────────────────────

def create_arg_yamls(exp_dict):

    NON_YAML_ARGS = ["demographies", "mating", "experiment", "custom_sims", "end_chr", "post_processing"]

    shared_args = {i: j for i, j in exp_dict.items() if i not in NON_YAML_ARGS}
    shared_args["base_dir"] = exp_dict["base_dir"]
    shared_args["post_process"] = exp_dict.get("post_processing", None)

    yaml_dir = Path(exp_dict["base_dir"]) / "yaml_files"
    all_end_chrs = list(exp_dict["end_chr"].keys())

    yaml_paths = []  # collect for status / commands

    # ── demographies x mating x end_chr ──────────────────────────────────────
    for demo_name, demo in exp_dict.get("demographies", {}).items():
        for mating_name, mating in exp_dict["mating"].items():
            for end_chr, end_chr_options in exp_dict["end_chr"].items():
                label = get_label(demo_name, mating_name, end_chr, all_end_chrs)
                resources = extract_resources(demo, mating, end_chr_options)
                exp_args = {
                    "custom_demo": demo,
                    "pedigree": mating,
                    "end_chr": end_chr,
                    "label": label,
                }
                yaml_path = yaml_dir / f"{label}.yaml"
                write_out_yaml(shared_args, exp_args, resources, yaml_path)
                yaml_paths.append(yaml_path)
                print(f"  wrote {yaml_path}")

    # ── custom_sims (no end_chr or mating iteration) ────────────────────────
    for custom_sim_name, custom_sim in exp_dict.get("custom_sims", {}).items():
        label = custom_sim_name
        resources = extract_resources(custom_sim)
        exp_args = {
            "custom_sim": custom_sim,
            "label": label,
        }
        yaml_path = yaml_dir / f"{label}.yaml"
        write_out_yaml(shared_args, exp_args, resources, yaml_path)
        yaml_paths.append(yaml_path)
        print(f"  wrote {yaml_path}")

    return yaml_paths


# ── Summary ───────────────────────────────────────────────────────────────────

def print_summary(exp_dict):
    """Print a human-readable overview of the experiment plan without creating anything."""

    exp_name = exp_dict["experiment"]
    demographies = exp_dict.get("demographies", {})
    custom_sims = exp_dict.get("custom_sims", {})
    mating = exp_dict.get("mating", {})
    end_chrs = list(exp_dict.get("end_chr", {}).keys())
    post_processing = exp_dict.get("post_processing", None)

    n_demo_runs = len(demographies) * len(mating) * len(end_chrs)
    n_custom_runs = len(custom_sims)
    n_total = n_demo_runs + n_custom_runs

    print(f"\nExperiment: {exp_name}")
    print("=" * 50)

    print(f"\nShared parameters:")
    print(f"  iter:     {exp_dict.get('iter')}")
    print(f"  samples:  {exp_dict.get('samples')}")
    print(f"  gb:       {exp_dict.get('gb')}")
    print(f"  sim_min:  {exp_dict.get('sim_min')}")
    print(f"  nthreads: {exp_dict.get('nthreads')}")

    print(f"\nAxes:")
    print(f"  end_chr:      {end_chrs}")
    print(f"  demographies: {list(demographies.keys())}")
    if custom_sims:
        print(f"  custom_sims:  {list(custom_sims.keys())} (no end_chr or mating iteration)")
    print(f"  mating:       {list(mating.keys())}")

    if demographies:
        print(f"\nDemography runs ({len(demographies)} demo x {len(end_chrs)} end_chr x {len(mating)} mating = {n_demo_runs} total):")
        for demo_name, demo in demographies.items():
            res = demo.get("resources", {})
            res_str = f"  [resources: {res}]" if res else ""
            print(f"  {demo_name}{res_str}")

    if custom_sims:
        print(f"\nCustom sim runs ({len(custom_sims)} total):")
        for sim_name, sim in custom_sims.items():
            res = sim.get("resources", {})
            res_str = f"  [resources: {res}]" if res else ""
            print(f"  {sim_name}{res_str}")

    print(f"\nMating models:")
    for mating_name in mating:
        print(f"  {mating_name}")

    print(f"\nPost-processing: {post_processing or 'none'}")
    print(f"\nTotal simulations: {n_total}")
    print()


# ── Status ────────────────────────────────────────────────────────────────────

def get_status(exp_dict, yaml_paths):
    """Check which simulations are complete, pending, or not started."""
    n_iter = exp_dict["iter"]
    base_dir = Path(exp_dict["base_dir"])

    rows = []
    for yaml_path in yaml_paths:
        args = yaml.safe_load(open(yaml_path))
        label = args["label"]

        # The output directory is base_dir/label (simulate.py uses label as dir name)
        out_dir = base_dir / label

        if not out_dir.exists():
            status = "not started"
            completed = 0
        else:
            completed = sum(
                1 for i in range(1, n_iter + 1)
                if (out_dir / f"iter{i}.ibd.gz").exists()
            )
            if completed == n_iter:
                status = "complete"
            elif completed > 0:
                status = f"in progress ({completed}/{n_iter})"
            else:
                status = "not started"

        rows.append((label, status, completed, n_iter))

    return rows


def print_status(exp_dict, yaml_paths):
    rows = get_status(exp_dict, yaml_paths)
    col_width = max(len(r[0]) for r in rows) + 2
    print(f"\n{'Simulation':<{col_width}} {'Status'}")
    print("-" * (col_width + 20))
    for label, status, completed, n_iter in rows:
        print(f"{label:<{col_width}} {status}")
    n_complete = sum(1 for _, status, _, _ in rows if status == "complete")
    print(f"\n{n_complete}/{len(rows)} simulations complete.\n")


# ── Commands ──────────────────────────────────────────────────────────────────

def print_commands(exp_dict, yaml_paths, pending_only=False, no_wait=False):
    """Print run.py commands for all (or only pending) simulations."""
    if pending_only:
        rows = get_status(exp_dict, yaml_paths)
        pending = {
            yaml_paths[i]
            for i, (_, status, _, _) in enumerate(rows)
            if status != "complete"
        }
        yaml_paths = [p for p in yaml_paths if p in pending]

    flag = " --no-wait" if no_wait else ""
    for yaml_path in yaml_paths:
        print(f"python run.py simulate {yaml_path}{flag}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        prog="experiment.py",
        description="IBD-sims experiment manager — plan and track batches of simulations.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # init
    p_init = sub.add_parser("init", help="Create experiment directory and generate YAML files")
    p_init.add_argument("yaml", help="Path to experiment meta-YAML")

    # status
    p_status = sub.add_parser("status", help="Print status of all simulations in an experiment")
    p_status.add_argument("yaml", help="Path to experiment meta-YAML")

    # commands
    p_cmd = sub.add_parser("commands", help="Print run.py commands for simulations")
    p_cmd.add_argument("yaml", help="Path to experiment meta-YAML")
    p_cmd.add_argument("--pending-only", action="store_true", default=False,
        help="Only print commands for simulations that are not yet complete")
    p_cmd.add_argument("--no-wait", action="store_true", default=False,
        help="Add --no-wait flag to generated commands")

    # describe
    p_desc = sub.add_parser("describe", help="Print a summary of the experiment plan (no files created)")
    p_desc.add_argument("yaml", help="Path to experiment meta-YAML")

    args = parser.parse_args()
    exp_dict = load_exp_yaml(args.yaml)

    if args.command == "init":
        print(f"Initialising experiment '{exp_dict['experiment']}'...")
        yaml_dir = setup_dirs(exp_dict)
        yaml_paths = create_arg_yamls(exp_dict)
        print(f"\nGenerated {len(yaml_paths)} YAML files in {yaml_dir}.")
        print(f"Run 'python experiment.py commands {args.yaml}' to get the commands to run.")

    elif args.command == "status":
        setup_dirs(exp_dict)
        yaml_dir = Path(exp_dict["base_dir"]) / "yaml_files"
        yaml_paths = sorted(yaml_dir.glob("*.yaml"))
        if not yaml_paths:
            print("No YAML files found. Run 'init' first.")
        else:
            print_status(exp_dict, yaml_paths)

    elif args.command == "commands":
        setup_dirs(exp_dict)
        yaml_dir = Path(exp_dict["base_dir"]) / "yaml_files"
        yaml_paths = sorted(yaml_dir.glob("*.yaml"))
        if not yaml_paths:
            print("No YAML files found. Run 'init' first.")
        else:
            print_commands(exp_dict, yaml_paths,
                           pending_only=args.pending_only,
                           no_wait=args.no_wait)


    elif args.command == "describe":
        print_summary(exp_dict)

if __name__ == "__main__":
    main()