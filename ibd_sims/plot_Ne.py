"""
plot_Ne.py — plot IBDNe and/or HapNe-IBD Ne estimates against the truth.

Usage
-----
python plot_Ne.py path/to/experiment_dir

Loads all results via load_experiment_results() and produces one Ne_plot.png
per demographic scenario found in the experiment directory.
"""

import argparse
import pickle
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

from simulations import DemographicSetup
from analyze_experiment import load_experiment_results
from demography import euro_bottleneck

QUEBEC_FOUNDING_GEN = 17  # BALSAC pedigree depth; Euro demography picks up here

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

def _is_quebec(yargs: dict) -> bool:
    """Return True if this demo uses the Quebec custom simulation."""
    custom_sim_path = (yargs.get("custom_sim") or {}).get("path") or ""
    return "quebec" in custom_sim_path.lower()


def get_truth(yargs: dict) -> pd.DataFrame | None:
    """Return a GEN/NE DataFrame representing the true demographic history.

    For Quebec: the BALSAC pedigree covers gens 0–QUEBEC_FOUNDING_GEN, so there
    is no analytic Ne truth for that window.  Uncoalesced lineages then continue
    through a European-like demography, so we return the euro_bottleneck
    trajectory starting at gen QUEBEC_FOUNDING_GEN.

    For all other custom_sim demos: no analytic truth available, return None.
    """
    gmax = yargs.get("gmax", 300)

    if _is_quebec(yargs):
        gen_arr = np.arange(QUEBEC_FOUNDING_GEN, gmax + 1)
        ne_arr = euro_bottleneck.debug().population_size_trajectory(gen_arr)
        return pd.DataFrame({"GEN": gen_arr, "NE": ne_arr[:, 0]})

    if (yargs.get("custom_sim") or {}).get("path"):
        return None  # other empirical/custom sims — no analytic truth

    try:
        demo = DemographicSetup.create(yargs)
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

OUTLIER_FACTOR = 100.0  # remove iters whose max NE exceeds this * median(max NE)


def _tool_key(label: str) -> str:
    """Infer tool name from a line label for palette selection."""
    low = label.lower()
    if "hapne-ibd" in low or "hapne_ibd" in low:
        return "hapne_ibd"
    if "hapne-ld" in low or "hapne_ld" in low:
        return "hapne_ld"
    return "ibdne"


def _filter_outlier_iters(
    dfs: list[pd.DataFrame],
    factor: float = OUTLIER_FACTOR,
) -> list[pd.DataFrame]:
    """
    Remove iterations whose max NE is more than `factor` times the median
    max NE across all iterations.  Protects against numerical blow-ups
    (e.g. IBDNe diverging on a single replicate) that would otherwise
    dominate the mean.
    """
    if len(dfs) <= 1:
        return dfs
    max_nes = np.array([df["NE"].max() for df in dfs])
    threshold = factor * np.median(max_nes)
    kept = [df for df, m in zip(dfs, max_nes) if m <= threshold]
    n_removed = len(dfs) - len(kept)
    if n_removed:
        print(f"    [outlier filter] Removed {n_removed}/{len(dfs)} iter(s) "
              f"(max NE > {factor:.0f}x median)")
    return kept if kept else dfs  # never drop everything


def _plot_panel(
    ax,
    panel_dict: dict,
    common_parts: list[str],
    palette_key: str,
    truth_df: pd.DataFrame | None,
    xlim: tuple,
    log_scale: bool,
    vlines: bool,
) -> None:
    """Plot one tool's data onto a single axes object."""
    line_styles = ["-", "--", "-.", ":"]
    palette = _TOOL_PALETTES.get(palette_key, _FALLBACK_PALETTE)

    for i, (label, dfs) in enumerate(panel_dict.items()):
        if not dfs:
            print(f"  No data for '{label}', skipping.")
            continue

        dfs = _filter_outlier_iters(dfs)

        color     = palette[i % len(palette)]
        linestyle = line_styles[i % len(line_styles)]

        # Strip parts shared across all lines (already in the panel title)
        legend_parts = [p for p in label.split("\n") if p not in common_parts]
        legend_label = "\n".join(legend_parts) if legend_parts else label

        # Align iteration lengths and stack
        ne_arrays = [df["NE"].values for df in dfs]
        min_len   = min(len(a) for a in ne_arrays)
        ne_arrays = [a[:min_len] for a in ne_arrays]
        generations = dfs[0]["GEN"].values[:min_len]
        ne_stack  = np.array(ne_arrays)   # (n_iter, n_gen)

        mean = np.mean(ne_stack, axis=0)
        ax.plot(generations, mean, label=legend_label, color=color,
                linewidth=2.5, linestyle=linestyle)

        if len(dfs) > BAND_THRESHOLD:
            p5  = np.percentile(ne_stack,  5, axis=0)
            p95 = np.percentile(ne_stack, 95, axis=0)
            ax.fill_between(generations, p5, p95, color=color, alpha=0.15)

    # Truth
    if truth_df is not None:
        ax.plot(truth_df["GEN"], truth_df["NE"],
                color="black", linestyle="--", label="Truth", linewidth=3.5)

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
    ax.set_xlabel("Generation")
    ax.grid(True, alpha=0.2)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left",
              borderaxespad=0., framealpha=1.0)


def plot_ne_estimates(
    data_dict: dict,
    truth_df: pd.DataFrame | None = None,
    figsize: tuple = (18, 7),
    log_scale: bool = True,
    xlim: tuple = (0, 50),
    vlines: bool = True,
) -> tuple:
    """
    Plot Ne estimates split into two panels (IBDNe | HapNe-IBD) with sharey=True.

    Parameters
    ----------
    data_dict : {label: [DataFrame, ...]}
        Each DataFrame must have GEN and NE columns.
    truth_df  : DataFrame with GEN and NE columns, or None.
    vlines    : draw log2-Ne vertical reference lines when truth is constant.
    """
    ibdne_dict = {l: d for l, d in data_dict.items() if _tool_key(l) == "ibdne"}
    hapne_dict = {l: d for l, d in data_dict.items() if _tool_key(l) == "hapne_ibd"}

    # Parts common to every label → use as suptitle, strip from legend entries
    all_parts = [label.split("\n") for label in data_dict] if data_dict else [[]]
    common = [p for p in all_parts[0]
              if all(p in parts for parts in all_parts[1:])]

    fig, (ax_left, ax_right) = plt.subplots(1, 2, sharey=True, figsize=figsize)
    fig.suptitle(" | ".join(common) if common else "", y=1.01)

    _plot_panel(ax_left,  ibdne_dict, common, "ibdne",     truth_df, xlim, log_scale, vlines)
    ax_left.set_title("IBDNe")
    ax_left.set_ylabel("Effective Population Size")

    _plot_panel(ax_right, hapne_dict, common, "hapne_ibd", truth_df, xlim, log_scale, vlines)
    ax_right.set_title("HapNe-IBD")

    plt.tight_layout()
    return fig, (ax_left, ax_right)


# ── Pickle helpers ────────────────────────────────────────────────────────────

_DISPLAY_TO_METHOD = {v: k for k, v in _METHOD_DISPLAY.items()}


def _parse_label(label: str) -> dict:
    """Parse a plot label string into a structured dict of metadata fields.

    Label format: "{demo}\\n{method_display} | rep={rep} | key=val | ..."
    """
    demo_part, rest = label.split("\n", 1)
    parts = [p.strip() for p in rest.split("|")]
    method_display = parts[0].strip()
    meta: dict = {
        "demo":   demo_part,
        "method": _DISPLAY_TO_METHOD.get(method_display, method_display),
    }
    for part in parts[1:]:
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        k, v = k.strip(), v.strip()
        # coerce booleans and numbers
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            for cast in (int, float):
                try:
                    v = cast(v)
                    break
                except ValueError:
                    pass
        meta[k] = v
    return meta


def _data_dict_to_records(data_dict: dict) -> list[dict]:
    """Convert the {label: [df, ...]} data_dict into a list of structured records.

    Each record has parsed metadata fields (demo, method, rep, any postprocess
    args) plus a 'dfs' key containing the list of per-iteration DataFrames.
    """
    records = []
    for label, dfs in data_dict.items():
        record = _parse_label(label)
        record["dfs"] = dfs
        records.append(record)
    return records


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

def plot(exp_dir: str, vlines: bool = True, save_pickle: bool = False) -> None:
    """
    Load all Ne estimates from an experiment directory and produce one
    Ne_plot.png per demographic scenario, saved into [exp_dir]/plots/.

    Parameters
    ----------
    exp_dir      : path to the experiment directory (passed to load_experiment_results)
    vlines       : draw log2-Ne vertical reference lines (only when truth is constant)
    save_pickle  : if True, write [exp_dir]/plots/plot_data.pkl with all data
                   needed to reproduce the plots in a Jupyter notebook
    """
    exp_dir = Path(exp_dir)
    plots_dir = exp_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

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

    # Collected for optional pickle output.
    # Structure:
    #   {
    #     "demos": {
    #       demo_name: {
    #         "data_dict": {label: [df, ...]},   # GEN / NE columns
    #         "truth_df":  pd.DataFrame | None,  # GEN / NE columns
    #       }
    #     },
    #     "plot_kwargs": {"log_scale": bool, "xlim": tuple, "vlines": bool},
    #   }
    pickle_payload: dict = {
        "demos": {},
        "plot_kwargs": {"log_scale": True, "xlim": (0, 50), "vlines": vlines},
    }

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
        out_path = plots_dir / f"{demo}_Ne_plot.png"
        plt.savefig(out_path, dpi=600)
        print(f"  Saved: {out_path}")
        plt.close(fig)

        if save_pickle:
            pickle_payload["demos"][demo] = {
                "records":  _data_dict_to_records(data_dict),
                "truth_df": truth_df,
            }

    if save_pickle:
        pkl_path = plots_dir / "Ne_data.pkl"
        with open(pkl_path, "wb") as fh:
            pickle.dump(pickle_payload, fh)
        print(f"\nPickle saved: {pkl_path}")


# ── RMSE table ────────────────────────────────────────────────────────────────

def _rmse_for_record(
    record: dict,
    truth_df: pd.DataFrame,
    gen_range: tuple[int, int],
    log_scale: bool,
) -> float | None:
    """
    Compute mean RMSE across all iters in a record against truth_df.

    RMSE is computed on log10(Ne) if log_scale=True (default, recommended
    because Ne spans orders of magnitude).  Truth is linearly interpolated
    at each iter's GEN values within gen_range.

    Returns None if truth_df is None or no GEN values fall in range.
    """
    if truth_df is None:
        return None

    g_lo, g_hi = gen_range
    truth_gen = truth_df["GEN"].values.astype(float)
    truth_ne  = truth_df["NE"].values.astype(float)

    rmses = []
    for df in record["dfs"]:
        gens = df["GEN"].values.astype(float)
        nes  = df["NE"].values.astype(float)

        mask = (gens >= g_lo) & (gens <= g_hi)
        if not mask.any():
            continue

        gens_in = gens[mask]
        nes_in  = nes[mask]

        # interpolate truth at the iter's GEN values
        truth_at = np.interp(gens_in, truth_gen, truth_ne)

        if log_scale:
            # guard against zeros / negatives before log
            ok = (nes_in > 0) & (truth_at > 0)
            if not ok.any():
                continue
            diff = np.log10(nes_in[ok]) - np.log10(truth_at[ok])
        else:
            diff = nes_in - truth_at

        rmses.append(np.sqrt(np.mean(diff ** 2)))

    return float(np.mean(rmses)) if rmses else None


_RECORD_META_KEYS = {"demo", "method", "rep", "dfs", "iter"}


def _select_records(
    records: list[dict],
    method: str,
    criteria: dict,
) -> list[dict]:
    """Return records matching method and all key=value pairs in criteria."""
    out = []
    for r in records:
        if r["method"] != method:
            continue
        if all(r.get(k) == v for k, v in criteria.items()):
            out.append(r)
    return out


_FILTER_ORDER = ["unfiltered", "random", "related", "unrelated"]


def _build_row_configs(filter_values: list[str]) -> list[dict]:
    """
    For each filter value produce two row configs:
      - filtersamples=False (IBDNe) / no filtersamples constraint (HapNe-IBD)
      - filtersamples=True  (IBDNe) / HapNe-IBD cells left empty

    'unfiltered' gets only a single row (filtersamples not applicable).
    """
    configs = []
    for fv in filter_values:
        configs.append({
            "label":       fv,
            "ibdne_crit":  {"filter": fv, "filtersamples": False},
            "hapne_crit":  {"filter": fv},
            "hapne_empty": False,
        })
        configs.append({
                "label":       r"\makecell[l]{" + fv + r"\\(filtersamples=True)}",
                "ibdne_crit":  {"filter": fv, "filtersamples": True},
                "hapne_crit":  {"filter": fv},
                "hapne_empty": True,
            })
    return configs


def compute_rmse_table(
    pickle_path: str | Path,
    gen_ranges: list[tuple[int, int]] | None = None,
    log_scale: bool = True,
    filter_values: list[str] | None = None,
) -> None:
    """
    Load Ne_data.pkl and print a LaTeX table of RMSE for IBDNe and HapNe-IBD.

    For each filter value two rows are produced per demographic scenario:
      - filtersamples=False for IBDNe; HapNe-IBD uses the same filter (no
        filtersamples distinction for HapNe-IBD)
      - filtersamples=True for IBDNe; HapNe-IBD cells are left empty

    RMSE is computed on log10(Ne) by default (log_scale=True).

    Parameters
    ----------
    pickle_path    : path to Ne_data.pkl produced by plot(..., save_pickle=True)
    gen_ranges     : list of (g_lo, g_hi) tuples; defaults to [(0, 50), (0, 10)]
    log_scale      : if True, compute RMSE on log10(Ne)
    filter_values  : filter values to include; defaults to all found in the data,
                     ordered by _FILTER_ORDER
    """
    if gen_ranges is None:
        gen_ranges = [(0, 50), (0, 10)]

    with open(pickle_path, "rb") as fh:
        payload = pickle.load(fh)

    demos_data = payload["demos"]
    scale_str  = r"$\log_{10} N_e$" if log_scale else r"$N_e$"

    # Detect filter values from data if not specified
    if filter_values is None:
        found = set()
        for demo_payload in demos_data.values():
            for r in demo_payload["records"]:
                fv = r.get("filter")
                if fv is not None:
                    found.add(fv)
        filter_values = [f for f in _FILTER_ORDER if f in found] + \
                        sorted(found - set(_FILTER_ORDER))

    ROW_CONFIGS = _build_row_configs(filter_values)

    methods = ["ibdne", "hapne_ibd"]
    method_crit_key = {"ibdne": "ibdne_crit", "hapne_ibd": "hapne_crit"}

    # ── collect results ────────────────────────────────────────────────────────
    # {demo: [ {(method, gen_range): rmse | None}, ... ]}  one dict per row config
    all_rows: list[tuple[str, str, dict, bool]] = []  # (demo, label, rmse_dict, hapne_empty)

    for demo, demo_payload in sorted(demos_data.items()):
        truth_df = demo_payload["truth_df"]
        records  = demo_payload["records"]

        for cfg in ROW_CONFIGS:
            rmse_dict: dict = {}
            for method in methods:
                crit = cfg[method_crit_key[method]]
                matching = _select_records(records, method, crit)
                if not matching:
                    for gr in gen_ranges:
                        rmse_dict[(method, gr)] = None
                    continue
                pooled = {"dfs": [df for r in matching for df in r["dfs"]]}
                for gr in gen_ranges:
                    rmse_dict[(method, gr)] = _rmse_for_record(
                        pooled, truth_df, gr, log_scale
                    )
            all_rows.append((demo, cfg["label"], rmse_dict, cfg["hapne_empty"]))

    # ── build LaTeX table ──────────────────────────────────────────────────────
    method_display = {"ibdne": "IBDNe", "hapne_ibd": "HapNe-IBD"}
    col_headers = [
        r"\makecell{" + method_display[m] + r"\\" + f"({g_lo}--{g_hi} gen)" + "}"
        for m in methods
        for (g_lo, g_hi) in gen_ranges
    ]
    n_cols = len(col_headers)
    col_spec = "ll" + "r" * n_cols

    lines = [
        r"\begin{table}[ht]",
        r"  \centering",
        f"  \\caption{{RMSE of {scale_str} Ne estimates vs.\ truth}}",
        r"  \label{tab:rmse}",
        f"  \\begin{{tabular}}{{{col_spec}}}",
        r"    \toprule",
        "    Demographic scenario & Filtering & " + " & ".join(col_headers) + r" \\",
        r"    \midrule",
    ]

    # group rows by demo; emit \midrule between demo groups, \hline between rows
    # within a group
    demos = list(dict.fromkeys(demo for demo, _, _, _ in all_rows))  # ordered unique
    for d_idx, demo in enumerate(demos):
        demo_rows = [(lbl, rd, he) for (dm, lbl, rd, he) in all_rows if dm == demo]
        for r_idx, (label, rmse_dict, hapne_empty) in enumerate(demo_rows):
            demo_cell = demo if r_idx == 0 else ""
            cells = []
            for m in methods:
                for gr in gen_ranges:
                    if hapne_empty and m == "hapne_ibd":
                        cells.append("")
                    else:
                        val = rmse_dict.get((m, gr))
                        cells.append(f"{val:.3f}" if val is not None else "---")
            lines.append(
                "    " + demo_cell + " & " + label + " & "
                + " & ".join(cells) + r" \\"
            )
            # \hline between rows within a demo group (not after the last row)
            if r_idx < len(demo_rows) - 1:
                lines.append(r"    \hline")

        # \midrule between demo groups (not after the last demo)
        if d_idx < len(demos) - 1:
            lines.append(r"    \midrule")

    lines += [
        r"    \bottomrule",
        r"  \end{tabular}",
        r"\end{table}",
    ]

    print("\n".join(lines))


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
    parser.add_argument("--save-pickle", action="store_true", default=False,
                        help="Write [exp_dir]/plots/plot_data.pkl with all data "
                             "needed to reproduce plots in a Jupyter notebook")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    plot(exp_dir=args.exp_dir, vlines=not args.no_vlines,
         save_pickle=args.save_pickle)