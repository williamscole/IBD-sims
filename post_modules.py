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
from filter_ibd import filter_ibd
from purple import readin_ibd


class PostProcessTest(PostProcessor):

    sub_config_key = "test"

    def _single_iter(self, iter_n):
        time.sleep(60)
        with open(Path(self.out_dir) / f"iter{iter_n}.txt", "w") as f:
            f.write(time.ctime())

    def execute(self):
        self._execute_helper()

        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop()


class PostProcessPurple(PostProcessor):
    sub_config_key = "purple_nodes"

    def execute(self):
        self._execute_helper()

        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop()

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

    def execute(self):
        self._execute_helper()

        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop()

    def _single_iter(self, iter_n):
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
                cfg.ibd_filter,
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

class PostProcessHapNeLD:

    sub_config_key = "hapne_ld"

    def _tmp_map(self, input_map):
        pass

    def _single_iter(self, iter_n):

        prefix = f"{self.path}/iter{iter_n}"

        map_file = hapne_tmp_map(f"{prefix}.map")

class PostProcessHapNeLD(PostProcessor):
    sub_config_key = "hapne_ld"
    resource_fields = ["local", "workers", "mem_gb", "time_min"]

    def execute(self):
        self._execute_helper()

        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop()

    def _tmp_map(self, input_map: str) -> tuple[str, str]:
        return hapne_tmp_map(input_map)

    def _single_iter(self, iter_n: int):
        from run_hapne import run_hapne_ld

        cfg = self._get_sub_config()
        prefix = f"{self.path}/iter{iter_n}"
        iter_out_dir = os.path.join(self.out_dir, f"iter{iter_n}")
        os.makedirs(iter_out_dir, exist_ok=True)

        run_hapne_ld(
            vcf_file=prefix,
            input_map=f"{prefix}.map",
            output_folder=iter_out_dir,
            population_name=getattr(cfg, "filter", "unfiltered"),
            end_chr=self.config.end_chr,
            workers=self._get_resource("workers")
        )


class PostProcessHapNeIBD(PostProcessor):
    sub_config_key = "hapne_ibd"
    resource_fields = ["local", "workers", "mem_gb", "time_min"]

    def execute(self):
        self._execute_helper()

        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop()

    def _single_iter(self, iter_n: int):
        from run_hapne import run_hapne_ibd

        cfg = self._get_sub_config()
        prefix = f"{self.path}/iter{iter_n}"
        iter_out_dir = os.path.join(self.out_dir, f"iter{iter_n}")
        os.makedirs(iter_out_dir, exist_ok=True)

        run_hapne_ibd(
            ibd_file=f"{prefix}.ibd.gz",
            input_map=f"{prefix}.map",
            output_folder=iter_out_dir,
            population_name=getattr(cfg, "filter", "unfiltered"),
            nb_samples=self.config.samples,
            end_chr=self.config.end_chr,
        )