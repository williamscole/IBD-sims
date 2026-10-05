"""
count_outlier_iters.py — tally iterations dropped by plot_Ne's outlier filter.

Usage
-----
python count_outlier_iters.py path/to/experiment_dir
python count_outlier_iters.py path/to/experiment_dir --csv outliers.csv
python count_outlier_iters.py path/to/experiment_dir --factor 50

plot_Ne.plot() drops, per plotted line, every iteration whose max NE exceeds
OUTLIER_FACTOR x median(max NE across that line's iterations).  This script
loads the experiment exactly the way plot() does, runs the *same*
_filter_outlier_iters() function on each line, and reports how many iterations
were removed per demographic history (with a per-line breakdown).

A "line" is one (method, rep) combination within a demo, i.e. one legend entry
in the plots (e.g. "IBDNe | rep=001 | filter=random").

Note: ne_results.py (RMSE tables) does NOT apply this filter, so these counts
describe the plots only.
"""

import argparse
import contextlib
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Allow running from anywhere (plot_Ne imports its siblings as top-level modules).
sys.path.insert(0, str(Path(__file__).resolve().parent))

from analyze_experiment import load_experiment_results  # noqa: E402
from plot_Ne import (  # noqa: E402
    OUTLIER_FACTOR,
    _HAPNE_FIXED_COLS,
    _IBDNE_FIXED_COLS,
    _filter_outlier_iters,
    _make_label,
    _results_to_data_dict,
)


def _iter_ids_by_label(ibdne_df: pd.DataFrame, hapne_df: pd.DataFrame, demo: str) -> dict:
    """Map each plot label to the iteration numbers behind its dfs list.

    Mirrors the grouping in plot_Ne._results_to_data_dict (groupby rep, then
    groupby iter, both sorted), so position i in the dfs list corresponds to
    position i here.
    """
    out: dict[str, list[int]] = {}
    for method, df, fixed in (
        ("ibdne", ibdne_df, _IBDNE_FIXED_COLS),
        ("hapne_ibd", hapne_df, _HAPNE_FIXED_COLS),
    ):
        if df.empty:
            continue
        subset = df[df["demo"] == demo]
        for rep, rep_group in subset.groupby("rep"):
            label = _make_label(demo, method, rep, rep_group.iloc[0], fixed)
            out[label] = [int(i) for i, _ in rep_group.groupby("iter")]
    return out


def _short_label(label: str) -> str:
    """Drop the leading demo line from a plot label."""
    return label.split("\n", 1)[1] if "\n" in label else label


def tally(exp_dir: str, factor: float = OUTLIER_FACTOR) -> pd.DataFrame:
    """Return one row per (demo, line) with iteration counts and removed iters."""
    results = load_experiment_results(Path(exp_dir))
    ibdne_df, hapne_df = results["ibdne"], results["hapne_ibd"]

    demos: set = set()
    for df in (ibdne_df, hapne_df):
        if not df.empty:
            demos |= set(df["demo"].unique())

    rows = []
    for demo in sorted(demos):
        # Silence the per-line "[ibdne/rep=...] n iter(s)" chatter from plot_Ne.
        with contextlib.redirect_stdout(io.StringIO()):
            data_dict = _results_to_data_dict(ibdne_df, hapne_df, demo)
        iter_ids = _iter_ids_by_label(ibdne_df, hapne_df, demo)

        for label, dfs in data_dict.items():
            ids = iter_ids[label]
            assert len(ids) == len(dfs), f"iter id / df mismatch for {label!r}"

            with contextlib.redirect_stdout(io.StringIO()):
                kept = _filter_outlier_iters(dfs, factor)
            kept_ids = {id(d) for d in kept}
            removed_idx = [i for i, d in enumerate(dfs) if id(d) not in kept_ids]

            # Informational only (the decision itself comes from the real filter).
            max_nes = np.array([d["NE"].max() for d in dfs])
            median_max = float(np.median(max_nes)) if len(max_nes) else float("nan")

            rows.append({
                "demo": demo,
                "line": _short_label(label),
                "n_iters": len(dfs),
                "n_removed": len(removed_idx),
                "removed_iters": ",".join(str(ids[i]) for i in removed_idx),
                "median_max_ne": median_max,
                "removed_max_ne": ",".join(f"{max_nes[i]:.3g}" for i in removed_idx),
            })

    return pd.DataFrame(rows)


def _print_report(df: pd.DataFrame, factor: float) -> None:
    if df.empty:
        print("No results found.")
        return

    print(f"\nOutlier filter: drop iters with max NE > {factor:g} x median(max NE) "
          f"within each line\n")

    for demo, g in df.groupby("demo", sort=True):
        total, removed = int(g["n_iters"].sum()), int(g["n_removed"].sum())
        print(f"=== {demo}: {removed} of {total} iteration(s) removed "
              f"across {len(g)} line(s) ===")
        show = g[["line", "n_iters", "n_removed", "removed_iters", "removed_max_ne"]]
        print(show.to_string(index=False))
        print()

    summary = (
        df.groupby("demo", sort=True)[["n_iters", "n_removed"]].sum().astype(int)
    )
    print("--- Summary (iterations summed over all lines in each demo) ---")
    print(summary.to_string())
    print(f"\nTOTAL removed: {int(df['n_removed'].sum())} of {int(df['n_iters'].sum())}")


def _parse_args():
    p = argparse.ArgumentParser(
        description="Count iterations removed by plot_Ne's outlier filter, per demo.",
    )
    p.add_argument("exp_dir", help="Experiment directory (as passed to plot_Ne.py)")
    p.add_argument("--factor", type=float, default=OUTLIER_FACTOR,
                   help=f"Outlier factor (default: plot_Ne.OUTLIER_FACTOR = {OUTLIER_FACTOR:g})")
    p.add_argument("--csv", metavar="PATH", default=None,
                   help="Also write the per-line table to this CSV file")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    table = tally(args.exp_dir, factor=args.factor)
    _print_report(table, args.factor)
    if args.csv and not table.empty:
        table.to_csv(args.csv, index=False)
        print(f"\nWrote {args.csv}")
