"""
simulate.py — orchestrator for IBD-sims pipeline.

Usage:
    # Slurm
    python simulate.py yaml_files/arg1.yaml

    # Local (8 parallel workers)
    python simulate.py yaml_files/arg1.yaml --local --workers 8

    # Resume a previous run
    python simulate.py path/to/existing/run/
"""

import argparse
import math
import os
import re
import sys
import shutil
import subprocess
import tempfile
import time
import yaml
import submitit
import itertools as it
from pathlib import Path

from simulations import sim, base_seed, load_config
from concat_tmrca import concat_tmrca
from purple import readin_ibd
from filter_ibd import filter_ibd, write_samples
import plot as plot_module
from post_process import postprocess
from utils import apply_overrides


# Terminal Slurm states that indicate a job will not recover
_FAILED_STATES = {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED"}


def _sacct_state(job_id):
    """Query Slurm accounting for the final state of a job.

    Used as a fallback when submitit reports UNKNOWN — this happens for
    array tasks that completed and aged out of squeue before submitit's
    own sacct query finds them.
    """
    result = subprocess.run(
        ["sacct", "-j", str(job_id), "--format=State", "--noheader", "-P", "--parsable2"],
        capture_output=True, text=True,
    )
    # First line is the job-level state (skip .batch and .extern steps)
    lines = [l.strip() for l in result.stdout.strip().splitlines() if l.strip()]
    return lines[0] if lines else "UNKNOWN"


def _release_if_held(job_id):
    """Check squeue; if the job is held, release it and return True."""
    result = subprocess.run(
        ["squeue", "--job", str(job_id), "--format=%r", "--noheader"],
        capture_output=True, text=True,
    )
    reason = result.stdout.strip()
    if reason and "held" in reason.lower():
        subprocess.run(["scontrol", "release", str(job_id)], capture_output=True)
        print(f"Auto-released held job {job_id} (reason: {reason})")
        return True
    return False


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_args(path):
    return yaml.safe_load(open(f"{path}/args.yaml"))


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


def is_sim_complete(path, iter_n, chrom):
    """Check whether a simulation job already completed successfully."""
    ibd_file = f"{path}/iter{iter_n}_chr{chrom}.ibd.gz"
    tmrca_file = f"{path}/iter{iter_n}_chr{chrom}.tmrca.pkl"
    map_file = f"{path}/iter{iter_n}_chr{chrom}.map"
    return os.path.exists(ibd_file) and os.path.exists(tmrca_file) and os.path.exists(map_file)


def is_post_complete(path, iter_n):
    """Check whether post-processing for an iteration already completed."""
    return os.path.exists(f"{path}/iter{iter_n}.ibd.gz")


def wait_for_jobs(jobs, path=None, end_chr=None):
    """Poll jobs until all finish. Auto-releases held jobs. Return list of failed jobs."""
    pending = list(jobs)
    failed = []
    while pending:
        still_pending = []
        for job in pending:
            state = _sacct_state(job.job_id)
            if state == "COMPLETED":
                pass
            elif state in _FAILED_STATES:
                print(f"Job {job.job_id} terminated with state {state}")
                failed.append(job)
            else:
                _release_if_held(job.job_id)
                still_pending.append(job)
        pending = still_pending
        if pending:
            time.sleep(10)
    return failed


# ── Job functions ─────────────────────────────────────────────────────────────

def run_pedigree(path, iter_n):
    args = load_args(path)
    from simulations import DemographicSetup, WFSetup
    demography = DemographicSetup.create(args)
    args["iteration_seed"] = base_seed(path, iter_n)
    args["iter_n"] = iter_n
    WFSetup.create(args, path, demography)


def run_simulation(path, iter_n, chrom):
    """Phase 2 (per-chr mode): simulate one (iteration, chromosome) pair."""
    sim(path, iter_n, chrom)


def run_simulation_iter(path, iter_n):
    """Phase 2 (per-iter mode): simulate all chromosomes for one iteration sequentially,
    then concatenate per-chromosome outputs into iteration-level files.

    Each chromosome is run in a subprocess so that any sys.exit() calls inside
    sim() only terminate the subprocess rather than the whole job.
    """
    from multiprocessing import Process

    args = load_args(path)
    end_chr = args["end_chr"]
    sim_workers = args.get("sim_workers", 1)
    prefix = f"{path}/iter{iter_n}"

    chroms_to_run = [c for c in range(1, end_chr + 1) if not is_sim_complete(path, iter_n, c)]
    skipped = end_chr - len(chroms_to_run)
    print(f"[iter {iter_n}] starting: end_chr={end_chr}, sim_workers={sim_workers}, "
          f"to simulate={len(chroms_to_run)}, skipped={skipped}")

    for i in range(0, len(chroms_to_run), sim_workers):
        batch = chroms_to_run[i:i + sim_workers]
        procs = []
        for chrom in batch:
            print(f"[iter {iter_n}] chr {chrom}: simulating...")
            p = Process(target=run_simulation, args=(path, iter_n, chrom))
            p.start()
            procs.append((chrom, p))
        for chrom, p in procs:
            p.join()
            if p.exitcode != 0:
                raise RuntimeError(f"Simulation failed for iter {iter_n} chr {chrom} (exit code {p.exitcode})")
            print(f"[iter {iter_n}] chr {chrom}: done")

    print(f"[iter {iter_n}] concatenating outputs...")
    concat_ibd_files(prefix, end_chr)
    concat_map_files(prefix, end_chr)
    if not os.path.exists(f"{prefix}.tmrca.gz"):
        concat_tmrca(path, iter_n, end_chr)
    if args.get("keep_all_files", False):
        concat_vcf(prefix, end_chr)
        concat_ibd_files(prefix, end_chr, hbd=True)
    write_samples(f"{prefix}.ibd.gz", args["samples"])

    print(f"[iter {iter_n}] cleaning up per-chromosome files...")
    remove_ibd_chr_files(prefix, end_chr, hbd=False)
    remove_ibd_chr_files(prefix, end_chr, hbd=True)
    remove_vcf_chr_files(prefix, end_chr)
    remove_map_chr_files(prefix, end_chr)
    print(f"[iter {iter_n}] done")

def concat_ibd_files(prefix, end_chr, hbd: bool = False):
    ext = "hbd" if hbd else "ibd"
    ibd_file = f"{prefix}.{ext}"
    if not os.path.exists(f"{ibd_file}.gz"):
        with open(ibd_file, "w") as out:
            for chrom in range(1, end_chr + 1):
                chr_file = f"{prefix}_chr{chrom}.{ext}.gz"
                result = subprocess.run(["zcat", chr_file], capture_output=True, text=True)
                for line in result.stdout.splitlines():
                    parts = line.split()
                    parts[4] = str(chrom)
                    out.write(" ".join(parts) + "\n")
        subprocess.run(["gzip", ibd_file], check=True)

def remove_ibd_chr_files(prefix, end_chr, hbd: bool = False):
    ext = "hbd" if hbd else "ibd"
    for chrom in range(1, end_chr + 1):
        chr_file = f"{prefix}_chr{chrom}.{ext}.gz"
        if os.path.exists(chr_file):
            os.remove(chr_file)

def concat_vcf(prefix, end_chr):
    vcf_file = f"{prefix}.vcf.gz"
    chr_files = [f"{prefix}_chr{chrom}.vcf.gz" for chrom in range(1, end_chr + 1)]
    if not os.path.exists(vcf_file) and all(os.path.exists(f) for f in chr_files):
        subprocess.run(
            ["bcftools", "concat"] + chr_files + ["-Oz", "-o", vcf_file],
            check=True
        )

def remove_vcf_chr_files(prefix, end_chr):
    for chrom in range(1, end_chr + 1):
        chr_file = f"{prefix}_chr{chrom}.vcf.gz"
        if os.path.exists(chr_file):
            os.remove(chr_file)

def concat_map_files(prefix, end_chr):
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

def remove_map_chr_files(prefix, end_chr):
    for chrom in range(1, end_chr + 1):
        chr_map = f"{prefix}_chr{chrom}.map"
        if os.path.exists(chr_map):
            os.remove(chr_map)

def run_post_processing(path, iter_n):
    """Phase 3: concatenate outputs and run post-processing"""
    args = load_args(path)
    config = load_config()
    end_chr = args["end_chr"]
    prefix = f"{path}/iter{iter_n}"

    # Concat files
    concat_ibd_files(prefix, end_chr)
    concat_map_files(prefix, end_chr)

    tmrca_file = f"{prefix}.tmrca.gz"
    if not os.path.exists(tmrca_file):
        concat_tmrca(path, iter_n, end_chr)

    if args.get("keep_all_files", False):
        concat_vcf(prefix, end_chr)
        concat_ibd_files(prefix, end_chr, hbd=True)

    # Write out samples for filtering
    write_samples(f"{prefix}.ibd.gz", args["samples"])

    # Post process
    postprocess(args, path=path, n_iter=None, iter_n=iter_n)

    # Delete intermediate files
    remove_ibd_chr_files(prefix, end_chr, hbd=False)
    remove_ibd_chr_files(prefix, end_chr, hbd=True)
    remove_vcf_chr_files(prefix, end_chr)
    remove_map_chr_files(prefix, end_chr)

    print(f"Simulation complete: iter {iter_n}")

# ── Orchestrator ──────────────────────────────────────────────────────────────

def run(yaml_path, local, n_workers, overrides=None, wait=True, max_n_slurm_jobs=1000):
    
    # Set up output directory
    if os.path.isdir(yaml_path) and os.path.exists(f"{yaml_path}/args.yaml"):
        path = yaml_path
    else:
        raw_args = yaml.safe_load(open(yaml_path))
        if overrides:
            raw_args = apply_overrides(raw_args, overrides)
        path = make_output_dir(yaml_path, raw_args)
        with open(f"{path}/args.yaml", "w") as f:
            yaml.dump(raw_args, f, default_flow_style=False)

    print(f"Output directory: {path}")

    args = load_args(path)
    n_iter = args["iter"]
    end_chr = args["end_chr"]
    pedigree_mode = args.get("pedigree") and args["pedigree"].get("pedigree_mode", False)
    sim_timeout = args["sim_min"]
    sim_workers = args.get("sim_workers", 1)
    pp_timeout = args.get("ibdne", {}).get("time_min") or args.get("time_min", 120)

    # Set up executor
    if local:
        executor = submitit.LocalExecutor(folder=f"{path}/slurm")
    else:
        executor = submitit.SlurmExecutor(folder=f"{path}/slurm")
        executor.update_parameters(use_srun=False, additional_parameters={"export": "ALL"})

    # ── Phase 1: Pedigree creation ────────────────────────────────────────────
    if pedigree_mode and not args["pedigree"].get("pedigree_file"):
        print("Submitting pedigree jobs...")
        if local:
            executor.update_parameters(timeout_min=20)
        else:
            executor.update_parameters(mem=4096, time=20, cpus_per_task=1)

        ped_iters = list(range(1, n_iter + 1))
        ped_jobs = executor.map_array(run_pedigree, [path] * len(ped_iters), ped_iters)

        failed = wait_for_jobs(ped_jobs)
        if failed:
            print(f"{len(failed)} pedigree jobs failed, aborting.")
            sys.exit(1)

    # ── Phase 2: Simulation ───────────────────────────────────────────────────
    per_chr_mode = n_iter < 3
    mode_str = "per-chromosome" if per_chr_mode else "per-iteration"
    print(f"Submitting simulation jobs ({mode_str} mode)...")

    if local:
        timeout = sim_timeout if per_chr_mode else sim_timeout * end_chr
        executor.update_parameters(timeout_min=timeout)
    else:
        if per_chr_mode:
            executor.update_parameters(
                mem=args["gb"] * 1024,
                time=sim_timeout,
                cpus_per_task=1,
                array_parallelism=100,
            )
        else:
            executor.update_parameters(
                mem=args["gb"] * 1024 * sim_workers,
                time=math.ceil(sim_timeout * end_chr / sim_workers),
                cpus_per_task=sim_workers,
            )

    sim_jobs = {}  # {iter_n: [job, ...]}

    if per_chr_mode:
        tasks = [
            (iter_n, chrom)
            for iter_n in range(1, n_iter + 1)
            for chrom in range(1, end_chr + 1)
            if not is_sim_complete(path, iter_n, chrom)
        ]
        if tasks:
            if len(tasks) <= max_n_slurm_jobs:
                jobs = executor.map_array(
                    run_simulation,
                    [path] * len(tasks),
                    [t[0] for t in tasks],
                    [t[1] for t in tasks],
                )
                for (iter_n, chrom), job in zip(tasks, jobs):
                    sim_jobs.setdefault(iter_n, []).append(job)
            else:
                if not wait:
                    print(f"Warning: --no-wait ignored because {len(tasks)} tasks exceeds "
                        f"max_n_slurm_jobs={max_n_slurm_jobs}. Batching required.")
                    wait = True

                n_iter_in_batch = max(1, max_n_slurm_jobs // 4 // end_chr)
                n_tasks_in_batch = n_iter_in_batch * end_chr
                batches = [tasks[i:i+n_tasks_in_batch] for i in range(0, len(tasks), n_tasks_in_batch)]

                print(f"Submitting {len(tasks)} tasks in {len(batches)} batches of "
                    f"{n_tasks_in_batch} ({n_iter_in_batch} iterations each)...")

                for i, batch in enumerate(batches):
                    print(f"Submitting batch {i+1}/{len(batches)}...")
                    jobs = executor.map_array(
                        run_simulation,
                        [path] * len(batch),
                        [t[0] for t in batch],
                        [t[1] for t in batch],
                    )
                    for (iter_n, chrom), job in zip(batch, jobs):
                        sim_jobs.setdefault(iter_n, []).append(job)

                    failed = wait_for_jobs(jobs)
                    if failed:
                        print(f"Warning: {len(failed)} jobs failed in batch {i+1}")
    else:
        iter_tasks = [
            iter_n for iter_n in range(1, n_iter + 1)
            if not all(is_sim_complete(path, iter_n, chrom) for chrom in range(1, end_chr + 1))
        ]
        if iter_tasks:
            jobs = executor.map_array(
                run_simulation_iter,
                [path] * len(iter_tasks),
                iter_tasks,
            )
            for iter_n, job in zip(iter_tasks, jobs):
                sim_jobs[iter_n] = [job]

    # for iter_n, jobs in sim_jobs.items():
    #     print(f"iter {iter_n}: {list(jobs.keys())} job ids: {[j.job_id for j in jobs.values()]}")

    time.sleep(15)

    # for iter_n, jobs in sim_jobs.items():
    #     for chrom, job in jobs.items():
    #         print(f"iter {iter_n} chr {chrom}: job_id={job.job_id} state={job.state}")

    # ── Early exit if --no-wait ───────────────────────────────────────────────
    if not wait:
        all_job_ids = [
            job.job_id
            for jobs in sim_jobs.values()
            for job in jobs
        ]
        print(f"Submitted {len(all_job_ids)} simulation jobs (--no-wait). Job IDs:")
        for job_id in all_job_ids:
            print(f"  {job_id}")
        print(f"To post-process later: python run.py postprocess {path}")
        return

    # ── Phase 3: Post-processing ──────────────────────────────────────────────
    print("Waiting for simulations and submitting post-processing...")
    if local:
        executor.update_parameters(timeout_min=pp_timeout)
    else:
        executor.update_parameters(
            mem=int(args["gb"] * 1.8 * 1024),
            time=pp_timeout,
            cpus_per_task=args.get("workers", args.get("nthreads", 1)),
            additional_parameters={}
        )

    post_jobs = {}
    for iter_n in range(1, n_iter + 1):
        print(f"iter {iter_n}: is_post_complete={is_post_complete(path, iter_n)}")
        if is_post_complete(path, iter_n):
            print(f"Skipping post-processing iter {iter_n} (already complete)")
            continue

        iter_jobs = sim_jobs.get(iter_n, [])
        if iter_jobs:
            failed = wait_for_jobs(iter_jobs)
            if failed:
                print(f"Warning: {len(failed)} simulation jobs failed for iter {iter_n}, skipping post-processing")
                continue

        post_jobs[iter_n] = executor.submit(run_post_processing, path, iter_n)

    # Wait for all post-processing to finish
    failed_post = wait_for_jobs(list(post_jobs.values()))
    if failed_post:
        print(f"Warning: {len(failed_post)} post-processing jobs failed")

    # ── Plotting ──────────────────────────────────────────────────────────────
    successful_post = [i for i in post_jobs if post_jobs[i] not in failed_post]
    if successful_post:
        print("Running plot...")
        plot_module.plot(path)

    print("Done.")


# ── Entry point ───────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description="IBD-sims pipeline orchestrator")
    parser.add_argument("yaml", help="Path to simulation YAML file, or existing run directory to resume")
    parser.add_argument("--local", action="store_true", help="Run locally instead of submitting to Slurm")
    parser.add_argument("--workers", type=int, default=os.cpu_count(),
                        help="Number of parallel workers for local execution (default: all available CPUs)")
    parser.add_argument("--no-wait", action="store_true", default=False,
                        help="Submit Slurm jobs and exit without waiting for them to finish")
    parser.add_argument("--max-jobs", type=int, default=1000,
                        help="Max number of slurm jobs to submit.")
    parser.add_argument(
        "--set", nargs="*", metavar="KEY=VALUE", default=None, action="append",
        help=(
            "Override YAML args without editing the file. "
            "Use dot notation for nested keys. "
            "Examples: --set iter=5  --set pedigree.mating=mono run_ibdne=false"
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.yaml, args.local, args.workers, overrides=args.set, wait=not args.no_wait, max_n_slurm_jobs=args.max_jobs)