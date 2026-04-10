"""
plot.py — plot IBDNe and/or HapNe Ne estimates against the truth.

Usage
-----
# Plot specific IBDNe and HapNe-IBD subdirectory runs:
python plot.py path/to/run --ibdne 001 003 --hapne_ibd 002

# Plot only HapNe-LD runs:
python plot.py path/to/run --hapne_ld 001 002

# Mix everything:
python plot.py path/to/run --ibdne 001 --hapne_ibd 001 --hapne_ld 001

The script reads args.yaml from each numbered subdirectory to build line labels,
and reads the top-level args.yaml to derive the truth Ne trajectory.
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

from simulations import DemographicSetup

BAND_THRESHOLD = 10   # show percentile band only when n_iters > this


# ── Label helpers ─────────────────────────────────────────────────────────────

def _build_base_label(yargs: dict) -> str:
    """Build a multiline base label from simulation args.yaml fields."""
    parts = []

    custom_demo = (yargs.get("custom_demo") or {}).get("object")
    custom_sim  = (yargs.get("custom_sim")  or {}).get("object")

    if custom_demo:
        parts.append(custom_demo)
    elif custom_sim:
        parts.append(custom_sim)

    n = yargs.get("samples")
    if n:
        parts.append(f"n={n}")

    ped = yargs.get("pedigree") or {}
    if ped.get("pedigree_mode"):
        mating = ped.get("mating", "di")
        parts.append(f"DTWF ({mating})")
    else:
        parts.append("coalescent")

    return "\n".join(parts)


def _ibdne_label(subdir_yargs: dict) -> str:
    """Label for an IBDNe line: base label + filter info.

    PostProcessor.dump_config() flattens the sub-config to the top level, so
    the subdir args.yaml has 'filter' and 'filtersamples' as top-level keys.
    """
    base = _build_base_label(subdir_yargs)

    # filter lives inside the ibdne: block in the subdir args.yaml
    ibdne_block = subdir_yargs.get("ibdne") or {}
    filter_val = ibdne_block.get("filter")
    filtersamples = ibdne_block.get("filtersamples") or subdir_yargs.get("filtersamples", False)

    if filter_val and str(filter_val) not in ("null", "", "none", "None", "unfiltered"):
        filter_str = str(filter_val)
    elif filtersamples:
        filter_str = "filtersamples"
    else:
        filter_str = "unfiltered"

    return f"{base}\nIBDNe | {filter_str}"


def _hapne_label(subdir_yargs: dict, tool: str) -> str:
    """Label for a HapNe line: base label + tool name + filter info."""
    base = _build_base_label(subdir_yargs)

    tool_cfg = subdir_yargs.get(tool) or {}
    filter_val = tool_cfg.get("filter") or subdir_yargs.get("filter") or "unfiltered"
    if filter_val in ("null", "", None):
        filter_val = "unfiltered"

    tool_display = {"hapne_ibd": "HapNe-IBD", "hapne_ld": "HapNe-LD"}.get(tool, tool)
    return f"{base}\n{tool_display} | {filter_val}"


# ── Data loaders ──────────────────────────────────────────────────────────────

def _load_ibdne_run(subdir: str) -> tuple[str, list[pd.DataFrame]]:
    """
    Load IBDNe results from a numbered subdirectory (e.g. ibdne/001/).

    Returns (label, list_of_dataframes). Each DataFrame has GEN and NE columns.
    """
    yargs_path = os.path.join(subdir, "args.yaml")
    if not os.path.exists(yargs_path):
        raise FileNotFoundError(f"No args.yaml found in {subdir}")

    yargs = yaml.safe_load(open(yargs_path))
    n_iter = yargs.get("iter") or yargs.get("n_iter")
    if n_iter is None:
        raise ValueError(f"Could not determine iter count from {yargs_path}")

    dfs = []
    for i in range(1, n_iter + 1):
        ne_file = os.path.join(subdir, f"iter{i}.ne")
        if os.path.exists(ne_file):
            dfs.append(pd.read_csv(ne_file, sep="\t"))
        else:
            print(f"  [IBDNe] Missing {ne_file}, skipping iteration {i}")

    label = _ibdne_label(yargs)
    return label, dfs


def _load_hapne_run(subdir: str, tool: str) -> tuple[str, list[pd.DataFrame]]:
    """
    Load HapNe-IBD or HapNe-LD results from a numbered subdirectory.

    HapNe output lives at: {subdir}/iter{i}/HapNe/hapne.csv
    Columns: TIME, Q0.025, Q0.25, Q0.5, Q0.75, Q0.975

    Returns (label, list_of_dataframes). Each DataFrame is normalised to
    have GEN (= TIME) and NE (= Q0.5) columns so it works with the
    shared plotting logic; the full quantile columns are preserved.
    """
    yargs_path = os.path.join(subdir, "args.yaml")
    if not os.path.exists(yargs_path):
        raise FileNotFoundError(f"No args.yaml found in {subdir}")

    yargs = yaml.safe_load(open(yargs_path))
    n_iter = yargs.get("iter") or yargs.get("n_iter")
    if n_iter is None:
        raise ValueError(f"Could not determine iter count from {yargs_path}")

    dfs = []
    for i in range(1, n_iter + 1):
        hapne_file = os.path.join(subdir, f"iter{i}", "HapNe", "hapne.csv")
        if os.path.exists(hapne_file):
            df = pd.read_csv(hapne_file)
            # Normalise column names for shared plotting code
            df = df.rename(columns={"TIME": "GEN", "Q0.5": "NE"})
            dfs.append(df)
        else:
            print(f"  [{tool}] Missing {hapne_file}, skipping iteration {i}")

    label = _hapne_label(yargs, tool)
    return label, dfs


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


# ── Main entry ────────────────────────────────────────────────────────────────

def plot(
    path: str,
    ibdne: list[str] | None = None,
    hapne_ibd: list[str] | None = None,
    hapne_ld: list[str] | None = None,
    vlines: bool = True,
) -> None:
    """
    Build and save the Ne plot.

    Parameters
    ----------
    path      : top-level run directory (contains args.yaml)
    ibdne     : list of subdirectory numbers to plot from ibdne/, e.g. ["001", "003"]
    hapne_ibd : list of subdirectory numbers to plot from hapne_ibd/
    hapne_ld  : list of subdirectory numbers to plot from hapne_ld/
    vlines    : draw log2-Ne vertical reference lines (only when truth is constant)
    """
    top_yargs_path = os.path.join(path, "args.yaml")
    if not os.path.exists(top_yargs_path):
        print(f"No args.yaml found at {path}")
        return

    top_yargs = yaml.safe_load(open(top_yargs_path))
    truth_df = get_truth(top_yargs)

    data_dict: dict[str, list[pd.DataFrame]] = {}

    tool_specs = [
        ("ibdne",     ibdne,     _load_ibdne_run),
        ("hapne_ibd", hapne_ibd, lambda s: _load_hapne_run(s, "hapne_ibd")),
        ("hapne_ld",  hapne_ld,  lambda s: _load_hapne_run(s, "hapne_ld")),
    ]

    for tool_name, run_ids, loader in tool_specs:
        if not run_ids:
            continue
        for run_id in run_ids:
            subdir = os.path.join(path, tool_name, run_id)
            if not os.path.isdir(subdir):
                print(f"  [{tool_name}] Subdirectory not found: {subdir}")
                continue
            try:
                label, dfs = loader(subdir)
            except Exception as e:
                print(f"  [{tool_name}/{run_id}] Failed to load: {e}")
                continue

            if not dfs:
                print(f"  [{tool_name}/{run_id}] No data found, skipping.")
                continue

            # Deduplicate labels (shouldn't happen, but just in case)
            unique_label = label
            suffix = 1
            while unique_label in data_dict:
                unique_label = f"{label} ({suffix})"
                suffix += 1

            data_dict[unique_label] = dfs
            print(f"  [{tool_name}/{run_id}] Loaded {len(dfs)} iteration(s): '{label}'")

    if not data_dict:
        print("No data to plot.")
        return

    fig, ax = plot_ne_estimates(data_dict, truth_df=truth_df, vlines=vlines)
    out_path = os.path.join(path, "Ne_plot.png")
    plt.savefig(out_path, dpi=600)
    print(f"Saved: {out_path}")
    plt.close(fig)


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args():
    parser = argparse.ArgumentParser(
        description="Plot IBDNe / HapNe Ne estimates.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python plot.py path/to/run --ibdne 001 003 --hapne_ibd 001
  python plot.py path/to/run --hapne_ld 001 002
        """,
    )
    parser.add_argument("path", help="Top-level run directory (contains args.yaml)")
    parser.add_argument("--ibdne",     nargs="+", metavar="RUN", default=None,
                        help="IBDNe subdirectory numbers to include (e.g. 001 003)")
    parser.add_argument("--hapne_ibd", nargs="+", metavar="RUN", default=None,
                        help="HapNe-IBD subdirectory numbers to include")
    parser.add_argument("--hapne_ld",  nargs="+", metavar="RUN", default=None,
                        help="HapNe-LD subdirectory numbers to include")
    parser.add_argument("--no-vlines", action="store_true", default=False,
                        help="Suppress log2-Ne vertical reference lines")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    plot(
        path=args.path,
        ibdne=args.ibdne,
        hapne_ibd=args.hapne_ibd,
        hapne_ld=args.hapne_ld,
        vlines=not args.no_vlines,
    )