"""
ne_results.py -- load and summarise saved Ne estimates (Ne_data.pkl).

Deliberately light: depends only on numpy, pandas and the standard library, so it
can be imported from a figures/analysis environment that does not have msprime,
submitit, etc.  Do NOT add imports of other ibd_sims modules here.
"""

import pickle
from pathlib import Path

import numpy as np
import pandas as pd


def load_ne_data(exp_dir, pickle_name="plots/Ne_data.pkl", halve_hapne=True):
    """
    Load [exp_dir]/plots/Ne_data.pkl (written by plot_Ne.py --save-pickle).

    HapNe reports Ne in haploid units; IBDNe and the simulated truth are diploid.
    With halve_hapne=True (default) every HapNe record's NE column is divided by 2.
    The payload gets hapne_halved=True afterwards, so loading a pickle regenerated
    with the fix applied at the source cannot halve twice.
    """
    with open(Path(exp_dir) / pickle_name, "rb") as fh:
        data = pickle.load(fh)

    if halve_hapne and not data.get("hapne_halved", False):
        for demo_payload in data["demos"].values():
            for rec in demo_payload["records"]:
                if "hapne" not in rec["method"]:
                    continue
                for df in rec["dfs"]:
                    assert isinstance(df, pd.DataFrame) and {"GEN", "NE"} <= set(df.columns)
                    df["NE"] = df["NE"] / 2
        data["hapne_halved"] = True

    return data

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

    g_lo = max(g_lo, truth_gen.min())   # e.g. Quebec truth only starts at g=17

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


_DEMO_DISPLAY = {
    "OOA2__DTWF_di":             "Out-of-Africa",
    "constant_Ne_10k__DTWF_di":  r"Constant $N_e = 10{,}000$",
    "constant_Ne_100k__DTWF_di": r"Constant $N_e = 100{,}000$",
    "quebec":                    "Quebec",
}


def compute_rmse_table(
    data: "dict | str | Path",
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

    payload = data if isinstance(data, dict) else load_ne_data(data)

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
        rf"  \caption{{RMSE of {scale_str} estimates vs.\ truth}}",
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
            demo_cell = _DEMO_DISPLAY.get(demo, demo) if r_idx == 0 else ""
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

