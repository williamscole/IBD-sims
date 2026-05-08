"""
plot_Ne.py — plot IBDNe and/or HapNe-IBD Ne estimates against the truth.

Usage
-----
python plot_Ne.py path/to/experiment_dir

Loads all results via load_experiment_results() and produces one Ne_plot.png
per demographic scenario found in the experiment directory.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

from simulations import DemographicSetup
from analyze_experiment import load_experiment_results

BAND_THRESHOLD = 10   # show percentile band only when n_iters > this


# ── Label helpers ─────────────────────────────────────────────────────────────

# Fixed columns in each method's DataFrame — everything else is a postprocess arg.
_IBDNE_FIXED_COLS  = {"demo", "rep", "iter", "gen", "ne", "lwr_95ci", "upr_95ci"}
_HAPNE_FIXED_COLS  = {"demo", "rep", "iter", "time",
                      "ne_q025", "ne_q25", "ne_q50", "ne_q75", "ne_q975"}

_METHOD_DISPLAY = {"ibdne": "IBDNe", "hapne_ibd": "HapNe-IBD"}


def _postprocess_suffix(row: pd.Series, fixed_cols: set) -> str:
    """Return a label fragment from any postprocessing arg columns present."""
    extra = {k: v for k, v in row.items()
             if k not in fixed_cols and pd.notna(v)}
    return " | ".join(f"{k}={v}" for k, v in extra.items()) if extra else "unfiltered"


def _make_label(demo: str, method: str, rep: str,
                row: pd.Series, fixed_cols: set) -> str:
    method_display = _METHOD_DISPLAY.get(method, method)
    pp = _postprocess_suffix(row, fixed_cols)
    return f"{demo}\n{method_display} | rep={rep} | {pp}"


# ── Truth ─────────────────────────────────────────────────────────────────────

def get_truth(yargs: dict) -> pd.DataFrame | None:
    """Return a GEN/NE DataFrame representing the true demographic history."""
    if (yargs.get("custom_sim") or {}).get("path"):
        return None  # empirical/custom sim — no analytic truth

    try:
        demo = DemographicSetup.create(yargs)
        gmax = yargs.get("gmax", 300)
        gen_arr = np.arange(0, gmax + 1)
        ne_arr = demo.debug().population_size_trajectory(gen_arr)
        return pd.DataFrame({"GEN": gen_arr, "NE": ne_arr[:, 0]})
    except Exception as e:
        print(f"  [truth] Could not compute truth Ne: {e}")
        return None


# ── Core plotting ─────────────────────────────────────────────────────────────

# Per-tool color palettes: each is a list of hex colours, light→dark.
# Lines within a tool cycle through these; dashes cycle independently.
_TOOL_PALETTES = {
    "ibdne":     ["#52b788", "#2d6a4f", "#1b4332", "#40916c", "#74c69d"],  # greens
    "hapne_ibd": ["#e76f51", "#f4a261", "#c1440e", "#e9c46a", "#b5501a"],  # oranges/reds
    "hapne_ld":  ["#5e60ce", "#4361ee", "#3a0ca3", "#7b2d8b", "#480ca8"],  # blue-purples
}
_FALLBACK_PALETTE = sns.color_palette("deep", 20)


def _tool_key(label: str) -> str:
    """Infer tool name from a line label for palette selection."""
    low = label.lower()
    if "hapne-ibd" in low:
        return "hapne_ibd"
    if "hapne-ld" in low:
        return "hapne_ld"
    return "ibdne"


def plot_ne_estimates(
    data_dict: dict,
    truth_df: pd.DataFrame | None = None,
    figsize: tuple = (12, 8),
    log_scale: bool = True,
    xlim: tuple = (0, 50),
    vlines: bool = True,
) -> tuple:
    """
    Plot Ne estimates for one or more methods/runs.

    Parameters
    ----------
    data_dict : {label: [DataFrame, ...]}
        Each DataFrame must have GEN and NE columns.
        HapNe DataFrames may additionally have Q0.025/Q0.975 columns.
    truth_df : DataFrame with GEN and NE columns, or None.
    vlines : if True, draw the log2-Ne and 1.77*log2-Ne vertical reference
             lines when the truth is a constant Ne.
    """
    fig, ax = plt.subplots(figsize=figsize)
    line_styles = ["-", "--", "-.", ":"]

    # Assign colours per tool, cycling within each tool's palette
    tool_counters: dict[str, int] = {}
    line_colors: list = []
    for label in data_dict:
        tk = _tool_key(label)
        idx = tool_counters.get(tk, 0)
        palette = _TOOL_PALETTES.get(tk, _FALLBACK_PALETTE)
        line_colors.append(palette[idx % len(palette)])
        tool_counters[tk] = idx + 1

    # Compute shared title from strings common to all labels
    all_parts = [label.split("\n") for label in data_dict]
    common = [p for p in all_parts[0] if all(p in parts for parts in all_parts[1:])]
    plot_title = " | ".join(common) if common else ""

    for i, ((label, dfs), color) in enumerate(zip(data_dict.items(), line_colors)):
        if not dfs:
            print(f"  No data for '{label}', skipping.")
            continue

        linestyle = line_styles[i % len(line_styles)]
        # Strip common parts from the per-line legend label
        legend_parts = [p for p in label.split("\n") if p not in common]
        legend_label = "\n".join(legend_parts) if legend_parts else label

        ne_arrays = [df["NE"].values for df in dfs]

        # Align lengths (HapNe runs may differ in GEN extent across iters)
        min_len = min(len(a) for a in ne_arrays)
        ne_arrays = [a[:min_len] for a in ne_arrays]
        generations = dfs[0]["GEN"].values[:min_len]

        ne_stack = np.array(ne_arrays)  # shape: (n_iter, n_gen)
        mean = np.mean(ne_stack, axis=0)

        ax.plot(generations, mean, label=legend_label, color=color,
                linewidth=2.5, linestyle=linestyle)

        n_iter = len(dfs)
        if n_iter > BAND_THRESHOLD:
            p5  = np.percentile(ne_stack, 5,  axis=0)
            p95 = np.percentile(ne_stack, 95, axis=0)
            ax.fill_between(generations, p5, p95, color=color, alpha=0.15)

    # Truth
    if truth_df is not None:
        ax.plot(truth_df["GEN"], truth_df["NE"],
                color="black", linestyle="--",
                label="Truth", linewidth=3.5)

        if vlines and np.std(truth_df["NE"]) == 0:
            ne0 = truth_df.iloc[0]["NE"]
            ax.axvline(x=np.log(ne0) / np.log(2),
                       linestyle="--", color="black", linewidth=1.5,
                       label=r"$\log_2 N_e$")
            ax.axvline(x=1.77 * np.log(ne0) / np.log(2),
                       linestyle=":", color="black", linewidth=2,
                       label=r"$1.77 \times \log_2 N_e$")

    ax.set_xlim(*xlim)
    if log_scale:
        ax.set_yscale("log")

    ax.set_title(plot_title, pad=20)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Effective Population Size")
    ax.grid(True, alpha=0.2)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left",
               borderaxespad=0., framealpha=1.0)
    plt.tight_layout()
    return fig, ax


# ── Data conversion ───────────────────────────────────────────────────────────

def _results_to_data_dict(
    ibdne_df: pd.DataFrame,
    hapne_df: pd.DataFrame,
    demo: str,
) -> dict[str, list[pd.DataFrame]]:
    """
    Convert the flat DataFrames from load_experiment_results into the
    {label: [per-iter DataFrame, ...]} format expected by plot_ne_estimates.

    Each (method, rep) combination becomes one label; iters become the list.
    Columns are normalised to GEN and NE for the shared plotting logic.
    """
    data_dict: dict[str, list[pd.DataFrame]] = {}

    # IBDNe
    if not ibdne_df.empty:
        subset = ibdne_df[ibdne_df["demo"] == demo]
        for rep, rep_group in subset.groupby("rep"):
            first_row = rep_group.iloc[0]
            label = _make_label(demo, "ibdne", rep, first_row, _IBDNE_FIXED_COLS)
            dfs = [
                iter_df[["gen", "ne"]]
                    .rename(columns={"gen": "GEN", "ne": "NE"})
                    .reset_index(drop=True)
                for _, iter_df in rep_group.groupby("iter")
            ]
            data_dict[label] = dfs
            print(f"  [ibdne/rep={rep}] {len(dfs)} iter(s) for demo '{demo}'")

    # HapNe-IBD
    if not hapne_df.empty:
        subset = hapne_df[hapne_df["demo"] == demo]
        for rep, rep_group in subset.groupby("rep"):
            first_row = rep_group.iloc[0]
            label = _make_label(demo, "hapne_ibd", rep, first_row, _HAPNE_FIXED_COLS)
            dfs = [
                iter_df[["time", "ne_q50"]]
                    .rename(columns={"time": "GEN", "ne_q50": "NE"})
                    .reset_index(drop=True)
                for _, iter_df in rep_group.groupby("iter")
            ]
            data_dict[label] = dfs
            print(f"  [hapne_ibd/rep={rep}] {len(dfs)} iter(s) for demo '{demo}'")

    return data_dict


# ── Main entry ────────────────────────────────────────────────────────────────

def plot(exp_dir: str, vlines: bool = True) -> None:
    """
    Load all Ne estimates from an experiment directory and produce one
    Ne_plot.png per demographic scenario.

    Parameters
    ----------
    exp_dir : path to the experiment directory (passed to load_experiment_results)
    vlines  : draw log2-Ne vertical reference lines (only when truth is constant)
    """
    exp_dir = Path(exp_dir)
    results = load_experiment_results(exp_dir)
    ibdne_df  = results["ibdne"]
    hapne_df  = results["hapne_ibd"]

    demos = set()
    if not ibdne_df.empty:
        demos |= set(ibdne_df["demo"].unique())
    if not hapne_df.empty:
        demos |= set(hapne_df["demo"].unique())

    if not demos:
        print("No results found.")
        return

    for demo in sorted(demos):
        print(f"\nPlotting demo: {demo}")
        data_dict = _results_to_data_dict(ibdne_df, hapne_df, demo)

        if not data_dict:
            print(f"  No data for '{demo}', skipping.")
            continue

        # Load truth from the demo's args.yaml
        demo_args_path = exp_dir / demo / "args.yaml"
        truth_df = None
        if demo_args_path.exists():
            yargs = yaml.safe_load(open(demo_args_path))
            truth_df = get_truth(yargs)
        else:
            print(f"  [truth] No args.yaml found at {demo_args_path}")

        fig, ax = plot_ne_estimates(data_dict, truth_df=truth_df, vlines=vlines)
        out_path = exp_dir / f"{demo}_Ne_plot.png"
        plt.savefig(out_path, dpi=600)
        print(f"  Saved: {out_path}")
        plt.close(fig)


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args():
    parser = argparse.ArgumentParser(
        description="Plot IBDNe / HapNe-IBD Ne estimates for an experiment.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  python plot_Ne.py path/to/experiment_dir",
    )
    parser.add_argument("exp_dir", help="Experiment directory (passed to load_experiment_results)")
    parser.add_argument("--no-vlines", action="store_true", default=False,
                        help="Suppress log2-Ne vertical reference lines")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    plot(exp_dir=args.exp_dir, vlines=not args.no_vlines)