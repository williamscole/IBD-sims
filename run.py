#!/usr/bin/env python
"""
run.py — single entry point for the IBD-sims pipeline.

Usage:
    python run.py simulate yaml_files/arg1.yaml --local --workers 8
    python run.py postprocess path/to/run/ --set ibdne.mincm=3
    python run.py plot path/to/run/ --ibdne 001 003 --hapne_ibd 001
"""

import argparse
import os
import sys
from pathlib import Path

# Add ibd_sims/ to the import path so internal imports resolve unchanged
sys.path.insert(0, str(Path(__file__).parent / "ibd_sims"))


def cmd_simulate(args):
    from simulate import run
    run(args.yaml, args.local, args.workers, overrides=args.set)


def cmd_postprocess(args):
    from post_process import AnalysisConfig, postprocess

    input_path = Path(args.input)
    if input_path.is_dir() and (input_path / "args.yaml").exists():
        yaml_file = input_path / "args.yaml"
    elif input_path.is_file():
        yaml_file = input_path
    else:
        sys.exit(f"Error: input not found: {args.input}")

    config = AnalysisConfig.from_yaml(yaml_file, overrides=args.set)
    postprocess(config, n_iter=config.n_iter, path=config.path, wait=not args.no_wait)


def cmd_plot(args):
    from plot_Ne import plot
    plot(
        path=args.path,
        ibdne=args.ibdne,
        hapne_ibd=args.hapne_ibd,
        hapne_ld=args.hapne_ld,
        vlines=not args.no_vlines,
    )


def main():
    parser = argparse.ArgumentParser(
        prog="run.py",
        description="IBD-sims pipeline — simulate, post-process, and plot.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ── simulate ──────────────────────────────────────────────────────────
    p_sim = sub.add_parser("simulate", help="Run simulations")
    p_sim.add_argument("yaml",
        help="Path to simulation YAML file, or existing run directory to resume")
    p_sim.add_argument("--local", action="store_true",
        help="Run locally instead of submitting to Slurm")
    p_sim.add_argument("--workers", type=int, default=os.cpu_count(),
        help="Number of parallel workers for local execution (default: all CPUs)")
    p_sim.add_argument("--set", nargs="*", metavar="KEY=VALUE",
        default=None, action="append",
        help="Override YAML parameters (e.g. --set iter=5 pedigree.mating=mono)")
    p_sim.set_defaults(func=cmd_simulate)

    # ── postprocess ───────────────────────────────────────────────────────
    p_pp = sub.add_parser("postprocess", help="Run post-processing analyses")
    p_pp.add_argument("input",
        help="Run directory (containing args.yaml) or path to a YAML config")
    p_pp.add_argument("--no-wait", action="store_true", default=False,
        help="Submit Slurm jobs and exit without waiting")
    p_pp.add_argument("--set", nargs="*", metavar="KEY=VALUE",
        default=None, action="append",
        help="Override config values (e.g. --set ibdne.nboots=100 local=false)")
    p_pp.set_defaults(func=cmd_postprocess)

    # ── plot ──────────────────────────────────────────────────────────────
    p_plot = sub.add_parser("plot", help="Plot Ne estimates vs truth")
    p_plot.add_argument("path",
        help="Top-level run directory (contains args.yaml)")
    p_plot.add_argument("--ibdne", nargs="+", metavar="RUN", default=None,
        help="IBDNe subdirectory numbers (e.g. 001 003)")
    p_plot.add_argument("--hapne_ibd", nargs="+", metavar="RUN", default=None,
        help="HapNe-IBD subdirectory numbers")
    p_plot.add_argument("--hapne_ld", nargs="+", metavar="RUN", default=None,
        help="HapNe-LD subdirectory numbers")
    p_plot.add_argument("--no-vlines", action="store_true", default=False,
        help="Suppress log2-Ne vertical reference lines")
    p_plot.set_defaults(func=cmd_plot)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()