"""
run_pipeline.py — orchestrator for IBD-sims pipeline.

Usage:
    # Slurm
    python run_pipeline.py yaml_files/arg1.yaml

    # Local (8 parallel workers)
    python run_pipeline.py yaml_files/arg1.yaml --local --workers 8

    # Resume a previous run
    python run_pipeline.py path/to/existing/run/
"""

import argparse
import os
import re
import sys
import shutil
import subprocess
import gzip
import yaml
import submitit
import numpy as np
import itertools as it
from pathlib import Path
from concurrent.futures import wait as futures_wait, ALL_COMPLETED

from simulations import (
    WFSetup, Simulation, DemographicSetup,
    run_hapibd, add_tmrca, base_seed
)
from write_vcf import write_vcf
from concat_tmrca import concat_tmrca
from purple import readin_ibd
from filter_ibd import filter_ibd


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_config():
    config_path = Path(__file__).parent / "setup.yaml"
    return yaml.safe_load(open(config_path))


def load_args(path):
    return yaml.safe_load(open(f"{path}/args.yaml"))


def parse_slurm_time(time_str):
    """Convert a Slurm --time string to minutes. Accepts HH:MM:SS or MM:SS."""
    time_str = time_str.replace("--time=", "")
    parts = time_str.strip().split(":")
    if len(parts) == 3:
        h, m, s = parts
        return int(h) * 60 + int(m) + int(s) // 60
    elif len(parts) == 2:
        m, s = parts
        return int(m) + int(s) // 60
    else:
        raise ValueError(f"Unrecognised time format: {time_str}")


def make_output_dir(yaml_path, args):
    """Derive output directory name from label, handle clashes."""
    label = args.get("label", "ibdne_sim")
    # Flatten multiline label
    label = "_".join(line.strip() for line in label.strip().splitlines() if line.strip())
    # Clean up characters that are awkward in directory names
    label = re.sub(r"\s+", "-", label)
    label = re.sub(r"[()'\"]", "", label)

    base_dir = args.get("base_dir")
    if base_dir and base_dir != "null":
        os.makedirs(base_dir, exist_ok=True)
        label = os.path.join(base_dir, label)

    # Increment suffix to avoid clashes
    path = label
    counter = 1
    while os.path.exists(path):
        path = f"{label}_{counter:03d}"
        counter += 1

    os.makedirs(path)
    os.makedirs(f"{path}/slurm")
    os.makedirs(f"{path}/errors")
    shutil.copy(yaml_path, f"{path}/args.yaml")

    return path


def is_complete(path, iter_n, chrom):
    """Check whether a simulation job already completed successfully."""
    ibd_file = f"{path}/iter{iter_n}_chr{chrom}.ibd.gz"
    if not os.path.exists(ibd_file):
        return False
    # Must have at least 10 rows (same check as simulations.py)
    try:
        with gzip.open(ibd_file, "rt") as f:
            lines = [f.readline() for _ in range(10)]
        return len([l for l in lines if l.strip()]) == 10
    except Exception:
        return False


def is_post_complete(path, iter_n):
    """Check whether post-processing for an iteration already completed."""
    return os.path.exists(f"{path}/iter{iter_n}.ibd.gz")


# ── Job functions ─────────────────────────────────────────────────────────────

def run_pedigree(path, iter_n):
    """Phase 1: create WF pedigree for one iteration."""
    args = load_args(path)
    args["iter_n"] = iter_n
    args["iteration_seed"] = base_seed(path, iter_n)

    demography = DemographicSetup.create(args)
    WFSetup.create(args, f"{path}/iter{iter_n}", demography)
    print(f"Pedigree created: iter {iter_n}")


def run_simulation(path, iter_n, chrom):
    """Phase 2: simulate one (iteration, chromosome) pair."""
    args = load_args(path)
    iteration_seed = base_seed(path, int(iter_n))
    args["iter_n"] = iter_n
    args["iteration_seed"] = iteration_seed
    args["seed"] = iteration_seed + chrom  # unique per chrom
    args["chrom"] = chrom

    prefix = f"{path}/iter{iter_n}_chr{chrom}"

    ts, rate, _, seed = Simulation.create(args, chrom, f"{path}/iter{iter_n}")

    write_vcf(ts, prefix, chrom, rate, seed)
    run_hapibd(prefix, args["gb"])
    add_tmrca(prefix, ts, dump=False)

    # Cleanup if successful
    if os.path.exists(f"{prefix}.ibd.gz"):
        with gzip.open(f"{prefix}.ibd.gz", "rt") as f:
            lines = [f.readline() for _ in range(10)]
        if len([l for l in lines if l.strip()]) == 10:
            for ext in [".vcf.gz", ".hbd.gz"]:
                if os.path.exists(f"{prefix}{ext}"):
                    os.remove(f"{prefix}{ext}")
            print(f"Success: iter {iter_n}, chr {chrom}")
            return
    
    # Write error if unsuccessful
    os.makedirs(f"{path}/errors", exist_ok=True)
    with open(f"{path}/errors/iter{iter_n}_chr{chrom}.err", "w") as f:
        f.write(f"No .ibd.gz file found for iter {iter_n}, chromosome {chrom}\n")
    raise RuntimeError(f"Simulation failed: iter {iter_n}, chr {chrom}")


def run_post_processing(path, iter_n):
    """Phase 3: concatenate outputs and run IBDNe for one iteration."""
    args = load_args(path)
    config = load_config()
    end_chr = args["end_chr"]
    prefix = f"{path}/iter{iter_n}"

    # Purple node matrix
    purple_file = f"{prefix}.npy"
    if not os.path.exists(purple_file):
        readin_ibd(f"{prefix}_chr1.ibd.gz", n_chrom=end_chr, file_to_write=purple_file)

    # Concatenate TMRCAs
    tmrca_file = f"{prefix}.tmrca.gz"
    if not os.path.exists(tmrca_file):
        concat_tmrca(path, iter_n, end_chr)

    # Concatenate IBD files
    ibd_file = f"{prefix}.ibd"
    if not os.path.exists(f"{ibd_file}.gz"):
        with open(ibd_file, "w") as out:
            for chrom in range(1, end_chr + 1):
                chr_file = f"{prefix}_chr{chrom}.ibd.gz"
                result = subprocess.run(
                    ["zcat", chr_file],
                    capture_output=True, text=True
                )
                for line in result.stdout.splitlines():
                    parts = line.split()
                    parts[4] = str(chrom)
                    out.write(" ".join(parts) + "\n")
                os.remove(chr_file)
        subprocess.run(["gzip", ibd_file], check=True)

    # Concatenate map files
    map_file = f"{prefix}.map"
    if not os.path.exists(map_file):
        with open(map_file, "w") as out:
            for chrom in range(1, end_chr + 1):
                chr_map = f"{prefix}_chr{chrom}.map"
                with open(chr_map) as f:
                    for line in f:
                        parts = line.split()
                        parts[0] = str(chrom)
                        out.write(" ".join(parts) + "\n")
                os.remove(chr_map)

    # Run IBDNe
    if args.get("run_ibdne"):
        dir_name = args.get("dir_name", "ibdne")
        out_dir = f"{path}/{dir_name}"
        os.makedirs(out_dir, exist_ok=True)
        shutil.copy(f"{path}/args.yaml", f"{out_dir}/args.yaml")

        gb = int(args["gb"] * 1.8)

        filtered_ibd = filter_ibd(f"{prefix}.ibd.gz", args["samples"], args.get("ibd_filter"))

        ibdne_cmd = (
            f"java -jar -Xmx{gb}g {config['ibdne_jar']}"
            f" map={prefix}.map"
            f" out={out_dir}/iter{iter_n}"
            f" nthreads={args['nthreads']}"
            f" filtersamples={str(args['filtersamples']).lower()}"
            f" npairs={args['npairs']}"
            f" nits={args['nits']}"
            f" nboots={args['nboots']}"
            f" mincm={args['mincm']}"
            f" trimcm={args['trimcm']}"
            f" gmin={args['gmin']}"
            f" gmax={args['gmax']}"
        )

        proc = subprocess.run(
            ibdne_cmd,
            input=filtered_ibd,
            shell=True,
            capture_output=True,
            text=True
        )

        if proc.returncode != 0:
            raise RuntimeError(f"IBDNe failed for iter {iter_n}:\n{proc.stderr}")

        print(f"Post-processing complete: iter {iter_n}")


# ── Orchestrator ──────────────────────────────────────────────────────────────

def get_executor(args, path, local, n_workers):
    if local:
        executor = submitit.LocalExecutor(folder=f"{path}/slurm")
        executor.update_parameters(timeout_min=120, workers=n_workers)
    else:
        executor = submitit.AutoExecutor(folder=f"{path}/slurm")
        executor.update_parameters(
            slurm_partition="batch",
            slurm_account=None
        )
    return executor


def run(yaml_path, local, n_workers):
    # Load args
    raw_args = yaml.safe_load(open(yaml_path))

    # Set up output directory
    if os.path.isdir(yaml_path) and os.path.exists(f"{yaml_path}/args.yaml"):
        path = yaml_path  # resuming existing run
    else:
        path = make_output_dir(yaml_path, raw_args)

    print(f"Output directory: {path}")

    args = load_args(path)
    n_iter = args["iter"]
    end_chr = args["end_chr"]
    pedigree_mode = args["pedigree"]["pedigree_mode"]

    sim_timeout = parse_slurm_time(args["sim_time"])
    ibdne_timeout = parse_slurm_time(args["ibdne_time"])

    executor = get_executor(args, path, local, n_workers)

    # ── Phase 1: Pedigree creation ────────────────────────────────────────────
    ped_jobs = {}
    if pedigree_mode:
        print("Submitting pedigree jobs...")
        executor.update_parameters(
            **({} if local else {"slurm_mem_gb": 4, "slurm_time": 20}),
            **({"timeout_min": 20} if local else {}),
            cpus_per_task=1
        )
        for iter_n in range(1, n_iter + 1):
            ped_jobs[iter_n] = executor.submit(run_pedigree, path, iter_n)

    # ── Phase 2: Simulation ───────────────────────────────────────────────────
    print("Submitting simulation jobs...")
    if not local:
        executor.update_parameters(
            slurm_mem_gb=args["gb"],
            slurm_time=sim_timeout,
            cpus_per_task=1,
            slurm_array_parallelism=100
        )

    sim_jobs = {}
    for iter_n in range(1, n_iter + 1):
        sim_jobs[iter_n] = {}
        for chrom in range(1, end_chr + 1):
            if is_complete(path, iter_n, chrom):
                print(f"Skipping iter {iter_n} chr {chrom} (already complete)")
                continue

            if not local and pedigree_mode:
                # Pass pedigree job as dependency
                executor.update_parameters(
                    slurm_additional_parameters={"dependency": f"afterok:{ped_jobs[iter_n].job_id}"}
                )

            sim_jobs[iter_n][chrom] = executor.submit(run_simulation, path, iter_n, chrom)

    # ── Phase 3: Post-processing (submit as each iteration completes) ─────────
    print("Waiting for simulations and submitting post-processing...")
    if not local:
        executor.update_parameters(
            slurm_mem_gb=int(args["gb"] * 1.8),
            slurm_time=ibdne_timeout,
            cpus_per_task=args["nthreads"],
            slurm_additional_parameters={}
        )

    post_jobs = {}
    for iter_n in range(1, n_iter + 1):
        if is_post_complete(path, iter_n):
            print(f"Skipping post-processing iter {iter_n} (already complete)")
            continue

        # Wait for all simulation jobs for this iteration
        iter_futures = list(sim_jobs[iter_n].values())
        if iter_futures:
            futures_wait(iter_futures, return_when=ALL_COMPLETED)

        # Check for any failures
        failed = [f for f in iter_futures if f.exception() is not None]
        if failed:
            print(f"Warning: {len(failed)} simulation jobs failed for iter {iter_n}, skipping post-processing")
            continue

        post_jobs[iter_n] = executor.submit(run_post_processing, path, iter_n)

    # Wait for all post-processing to finish
    futures_wait(list(post_jobs.values()), return_when=ALL_COMPLETED)

    # ── Plotting ──────────────────────────────────────────────────────────────
    if any(f.exception() is None for f in post_jobs.values()):
        print("Running plot...")
        import plot
        plot.main(path)

    print("Done.")


# ── Entry point ───────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description="IBD-sims pipeline orchestrator")
    parser.add_argument("yaml", help="Path to simulation YAML file, or existing run directory to resume")
    parser.add_argument("--local", action="store_true", help="Run locally instead of submitting to Slurm")
    parser.add_argument("--workers", type=int, default=os.cpu_count(),
                        help="Number of parallel workers for local execution (default: all available CPUs)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.yaml, args.local, args.workers)