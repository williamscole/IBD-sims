import numpy as np
import os
import subprocess
import tempfile
import time
from pathlib import Path
import shutil
import pandas as pd
from configparser import ConfigParser

from post_process import PostProcessor
from simulations import load_config
from filter_ibd import filter_ibd, get_nodes, get_node_file_path
from purple import readin_ibd


class PostProcessTest(PostProcessor):

    sub_config_key = "test"

    def _single_iter(self, iter_n):
        time.sleep(60)
        with open(Path(self.out_dir) / f"iter{iter_n}.txt", "w") as f:
            f.write(time.ctime())

    def execute(self, wait=True):
        self._execute_helper()

        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop(wait=wait)


class PostProcessPurple(PostProcessor):
    sub_config_key = "purple_nodes"

    def execute(self, wait=True):
        self._execute_helper()

        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop(wait=wait)

    def _single_iter(self, iter_n):
        end_chr = self.config.end_chr
        prefix = f"{self.path}/iter{iter_n}"
        out_file = f"{self.out_dir}/iter{iter_n}.npy"

        if os.path.exists(out_file):
            print(f"Purple file already exists: {out_file}, skipping.")
            return

        readin_ibd(
            f"{prefix}_chr1.ibd.gz",
            n_chrom=end_chr,
            file_to_write=out_file,
        )


class PostProcessIBDNe(PostProcessor):
    sub_config_key = "ibdne"

    def is_iter_complete(self, iter_n: int) -> bool:
        out_file = f"{self.out_dir}/iter{iter_n}.ne"
        return os.path.exists(out_file) and sum(1 for _ in open(out_file)) > 10

    def execute(self, wait=True):
        self._execute_helper()

        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop(wait=wait)

    def _single_iter(self, iter_n):
        out_file = f"{self.out_dir}/iter{iter_n}.ne"
        if os.path.exists(out_file) and sum(1 for _ in open(out_file)) > 10:
            print(f"Already complete: {out_file}, skipping.")
            return

        cfg = self._get_sub_config()
        setup = load_config()
        prefix = f"{self.path}/iter{iter_n}"
        workers = self._get_resource("workers")
        mem_gb = self._get_resource("mem_gb")

        with tempfile.NamedTemporaryFile(suffix=".ibd", delete=False) as tmp:
            tmp_path = tmp.name

        try:
            filter_ibd(
                f"{prefix}.ibd.gz",
                self.config.samples,
                tmp_path,
                _get_filter(self),
            )

            ibdne_cmd = (
                f"cat {tmp_path} |"
                f" java -jar -Xmx{int(mem_gb * 0.8)}g {setup['ibdne_jar']}"
                f" map={prefix}.map"
                f" out={self.out_dir}/iter{iter_n}"
                f" nthreads={workers}"
                f" filtersamples={str(cfg.filtersamples).lower()}"
                f" npairs={cfg.npairs}"
                f" nits={cfg.nits}"
                f" nboots={cfg.nboots}"
                f" mincm={cfg.mincm}"
                f" trimcm={cfg.trimcm}"
                f" gmin={cfg.gmin}"
                f" gmax={cfg.gmax}"
            )

            proc = subprocess.run(ibdne_cmd, shell=True, capture_output=True, text=True)

            if proc.returncode != 0:
                raise RuntimeError(f"IBDNe failed for iter {iter_n}:\n{proc.stderr}")

        finally:
            os.remove(tmp_path)


def _needs_filtering(filtering):
    """Return True if the filtering value indicates actual sample subsetting."""
    return filtering is not None and filtering not in ("null", "", "none", "unfiltered")


def _get_filter(processor):
    """Resolve the filter value: sub-config overrides top-level.

    Checks processor's sub-config for 'filter' first; if not set,
    falls back to self.config.filter (the top-level default).
    """
    sub = processor._get_sub_config()
    value = getattr(sub, "filter", None)
    if value is not None:
        return value
    return getattr(processor.config, "filter", None)


class PostProcessHapNeLD(PostProcessor):
    sub_config_key = "hapne_ld"
    resource_fields = ["local", "workers", "mem_gb", "time_min"]

    def execute(self, wait=True):
        self._execute_helper()

        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop(wait=wait)

    def _tmp_map(self, input_map: str) -> tuple[str, str]:
        return hapne_tmp_map(input_map)

    def _single_iter(self, iter_n: int):
        from run_hapne import run_hapne_ld

        prefix = f"{self.path}/iter{iter_n}"
        iter_out_dir = os.path.join(self.out_dir, f"iter{iter_n}")
        # TODO: add completion check once HapNe-LD output format is known
        os.makedirs(iter_out_dir, exist_ok=True)

        filtering = _get_filter(self)

        # Resolve keep file: use the cached iter{i}_{label}.txt node file
        # (one ID per row), which is exactly the format HapNe expects.
        keep_file = None
        if _needs_filtering(filtering):
            ibd_path = f"{prefix}.ibd.gz"
            # Ensure the node file exists (computes + caches if needed)
            get_nodes(ibd_path, self.config.samples, filtering, getattr(self.config, "subsample_frac", 0.25))
            keep_file = get_node_file_path(ibd_path, filtering)
            print(f"[HapNe-LD] using keep file: {keep_file}")

        population_name = filtering if _needs_filtering(filtering) else "unfiltered"

        thin = getattr(self._get_sub_config(), "thin", None)

        run_hapne_ld(
            vcf_file=prefix,
            input_map=f"{prefix}.map",
            output_folder=iter_out_dir,
            population_name=population_name,
            end_chr=self.config.end_chr,
            workers=self._get_resource("workers"),
            keep_file=keep_file,
            thin=thin,
        )


class PostProcessHapNeIBD(PostProcessor):
    sub_config_key = "hapne_ibd"
    resource_fields = ["local", "workers", "mem_gb", "time_min"]

    def is_iter_complete(self, iter_n: int) -> bool:
        out_file = os.path.join(self.out_dir, f"iter{iter_n}", "HapNe", "hapne.csv")
        return os.path.exists(out_file) and sum(1 for _ in open(out_file)) > 10

    def execute(self, wait=True):
        self._execute_helper()

        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop(wait=wait)

    def _single_iter(self, iter_n: int):
        from run_hapne import run_hapne_ibd

        prefix = f"{self.path}/iter{iter_n}"
        iter_out_dir = os.path.join(self.out_dir, f"iter{iter_n}")
        out_file = os.path.join(iter_out_dir, "HapNe", "hapne.csv")
        if os.path.exists(out_file) and sum(1 for _ in open(out_file)) > 10:
            print(f"Already complete: {out_file}, skipping.")
            return

        os.makedirs(iter_out_dir, exist_ok=True)

        filtering = _get_filter(self)
        population_name = filtering if _needs_filtering(filtering) else "unfiltered"

        if _needs_filtering(filtering):
            # Filter IBD segments to a tmp file, then pass to HapNe-IBD
            tmp_dir = tempfile.mkdtemp(prefix="hapne_ibd_")
            tmp_ibd = os.path.join(tmp_dir, f"iter{iter_n}.ibd")

            try:
                filter_ibd(
                    f"{prefix}.ibd.gz",
                    self.config.samples,
                    tmp_ibd,
                    filtering,
                )

                if os.path.getsize(tmp_ibd) == 0:
                    print(f"Warning: no IBD segments remain after '{filtering}' filtering for iter{iter_n}, skipping.")
                    return

                # nb_samples for HapNe-IBD must reflect the filtered count
                nodes = get_nodes(f"{prefix}.ibd.gz", self.config.samples, filtering, getattr(self.config, "subsample_frac", 0.25))
                nb_samples = len(nodes)

                run_hapne_ibd(
                    ibd_file=tmp_ibd,
                    input_map=f"{prefix}.map",
                    output_folder=iter_out_dir,
                    population_name=population_name,
                    nb_samples=nb_samples,
                    end_chr=self.config.end_chr,
                )
            finally:
                shutil.rmtree(tmp_dir)
        else:
            run_hapne_ibd(
                ibd_file=f"{prefix}.ibd.gz",
                input_map=f"{prefix}.map",
                output_folder=iter_out_dir,
                population_name=population_name,
                nb_samples=self.config.samples,
                end_chr=self.config.end_chr,
            )


class PostProcessIBDSummary(PostProcessor):
    """Summarise IBD segments across iterations and filter modes.

    For each iteration and each configured filter mode, computes:
      - n_samples         : number of samples in the (sub)set
      - n_segments        : number of IBD segments surviving filtering
      - n_sharing_pairs   : number of pairs sharing at least one IBD segment
      - mean_pairwise_ibd_cM : mean total IBD (cM) across sharing pairs
      - total_ibd_cM      : total IBD (cM) across all sharing pairs

    Per-iteration results are written to iter{n}.tsv, then concatenated
    into a single ibd_summary.tsv once all iterations are complete.
    """

    sub_config_key = "ibd_summary"

    def is_iter_complete(self, iter_n: int) -> bool:
        return os.path.exists(os.path.join(self.out_dir, f"iter{iter_n}.tsv"))

    def execute(self, wait=True):
        self._execute_helper()

        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop(wait=wait)
            self._concatenate()

    def _single_iter(self, iter_n: int):
        out_file = os.path.join(self.out_dir, f"iter{iter_n}.tsv")
        if os.path.exists(out_file):
            print(f"Already complete: {out_file}, skipping.")
            return

        cfg = self._get_sub_config()
        filters = getattr(cfg, "filters", [None, "related", "unrelated"])
        ibd_path = f"{self.path}/iter{iter_n}.ibd.gz"

        # Load IBD file once; reuse for all filter modes
        ibd_df = pd.read_csv(ibd_path, sep=r"\s+", header=None)

        # Apply minimum segment length threshold (col 7 = length in cM)
        mincm = getattr(cfg, "mincm", 2)
        ibd_df = ibd_df[ibd_df[7] > mincm]

        rows = []
        for filtering in filters:
            if _needs_filtering(filtering):
                nodes = get_nodes(ibd_path, self.config.samples, filtering, getattr(self.config, "subsample_frac", 0.25))
                filtered = ibd_df[ibd_df[0].isin(nodes) & ibd_df[2].isin(nodes)]
                n_samples = len(nodes)
                filter_label = filtering
            else:
                filtered = ibd_df
                n_samples = self.config.samples
                filter_label = "unfiltered"

            n_segments = len(filtered)

            if n_segments > 0:
                pair_ibd = filtered.groupby([0, 2])[7].sum()
                n_sharing_pairs = len(pair_ibd)
                mean_pairwise_ibd_cM = pair_ibd.mean()
                total_ibd_cM = pair_ibd.sum()
            else:
                n_sharing_pairs = 0
                mean_pairwise_ibd_cM = float("nan")
                total_ibd_cM = 0.0

            rows.append({
                "iter": iter_n,
                "filter": filter_label,
                "n_samples": n_samples,
                "n_segments": n_segments,
                "n_sharing_pairs": n_sharing_pairs,
                "mean_pairwise_ibd_cM": mean_pairwise_ibd_cM,
                "total_ibd_cM": total_ibd_cM,
            })

        pd.DataFrame(rows).to_csv(out_file, sep="\t", index=False)
        print(f"Wrote {out_file}")

    def _concatenate(self):
        """Concatenate all per-iter TSVs into a single ibd_summary.tsv."""
        parts = []
        for iter_n in range(1, self.n_iter + 1):
            f = os.path.join(self.out_dir, f"iter{iter_n}.tsv")
            if os.path.exists(f):
                parts.append(pd.read_csv(f, sep="\t"))

        if parts:
            out = os.path.join(self.out_dir, "ibd_summary.tsv")
            pd.concat(parts, ignore_index=True).to_csv(out, sep="\t", index=False)
            print(f"Wrote {out}")