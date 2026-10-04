#!/usr/bin/env python
"""
Re-create Fig 2 (pair_tmrca): for IBD segments grouped by their most likely
age g, show which pair-TMRCA classes contributed them.

Inputs (from an IBD-sims run directory):
    iter{N}.ibd.gz    id1 hap1 id2 hap2 chr start end cM   (no header, whitespace)
    iter{N}.tmrca.gz  tab-separated, header; row i  <->  row i of the ibd file.
                      'tmrca' / 'proportion' are stringified lists: the time of
                      the true MRCA node(s) overlapping the detected segment and
                      the fraction of the segment each one covers.

Definitions used here (the original code was lost, so these are reconstructions
-- the validation printout compares them with numbers quoted in the paper):
    segment TMRCA   time of the MRCA node covering most of the segment
    pair TMRCA      minimum segment TMRCA over ALL segments of that pair
                    (unordered pair of individuals, all chromosomes), rounded
                    to whole generations
    most likely age argmax_g of Eq. 1 (Browning & Browning 2015) over g = 1..100,
                    given segment length l (cM) and the sampled population's
                    N(g).  gamma_0 is a constant in g, so it drops out of argmax.

Iterations: with no --iter, every iter{N}.ibd.gz that has a matching
iter{N}.tmrca.gz in the run directory is used.  Everything (pair TMRCA, ages) is
computed WITHIN each iteration -- individual names like tsk_807 are reused across
iterations -- and only the segment counts are pooled.  The pooled proportions are
plotted; the per-iteration table and across-iteration SD are written to CSV.

Two-panel figure (--two-panel):
    A  the reviewer-suggested view: P(segment TMRCA bin | segment length), TMRCA binned as
       1, 2, 3-5, 6-10, 11-20, 21-50, >50 generations.  By default the simulated segments are
       binned by length (--len-bins) and each segment's own TMRCA is read from tmrca.gz.
       Those IBD files start at ~2 cM, so there is no 1 cM bar unless you re-detect lower.
       With --a-analytical the same bars come from Eq. 1 and the N(g) trajectory alone: Eq. 1
       is integrated over each length bin, weighting lengths by the model's own length
       density (closed form), normalised over g = 1..--a-max-g.  No simulated segments are
       used, so it does not depend on the input data.  Both versions are always printed,
       with their difference, as a check of Eq. 1 against the simulation.
    B  the Fig 2 breakdown (--panel-b full|grouped).  Segments whose pair TMRCA is > 10 are
       shown as a grey ">10" class by default (--over10 drop|clip for the alternatives).

Usage:
    python make_fig2.py --run main_experiment/OOA2__DTWF_mono \
        --demo-file /path/to/IBD-sims/ibd_sims/demography.py --demo-object ooa2 \
        --pop pop_0 --out fig2_new                 # all iterations
    ... --iter 1                                   # or just one / a few: --iter 1 2 5
    ... --two-panel                                # also write fig2_new_twopanel.png
"""
import argparse
import importlib.util
import math
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MAX_G = 100        # ages searched when maximising Eq. 1
SHOW_G = 10        # x-axis: most likely ages 1..SHOW_G
FLOAT_RE = re.compile(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?")


# ----------------------------------------------------------------- demography
def load_ne_trajectory(demo_file, demo_object, pop_name, max_g=MAX_G):
    """N[g] for g = 1..max_g (index 0 is g = 1), as wf_pedigree.py computes it."""
    spec = importlib.util.spec_from_file_location("demo_module", demo_file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    demography = getattr(mod, demo_object)
    dbg = demography.debug()
    names = [p.name for p in dbg.demography.populations]
    if pop_name not in names:
        sys.exit(f"Population {pop_name!r} not in demography; available: {names}")
    idx = names.index(pop_name)
    return dbg.population_size_trajectory(np.arange(1, max_g + 1))[:, idx].astype(float)


# ---------------------------------------------------------------------- Eq. 1
def _eq1_terms(N):
    """g = 1..G and the length-independent part of log Eq. 1 (everything except -l*g/50).

    P(g | l, N) = (1/gamma_0) (g/50)^2 exp(-l g / 50) prod_{g'=1}^{g-1} (1 - 1/(2N[g'])) / (2N[g])
    """
    g = np.arange(1, len(N) + 1, dtype=float)
    log_surv = np.concatenate([[0.0], np.cumsum(np.log1p(-1.0 / (2 * N)))[:-1]])   # prod over g' < g
    return g, 2 * np.log(g / 50.0) + log_surv - np.log(2 * N)


def most_likely_age(lengths_cm, N):
    """argmax_g P(TMRCA = g | l, N) for each length (vectorised over unique l).

    gamma_0 is constant in g, so it drops out of the argmax.
    """
    g, base = _eq1_terms(N)                                     # base: (G,)

    uniq, inv = np.unique(np.round(lengths_cm, 3), return_inverse=True)
    best = np.empty(len(uniq), dtype=int)
    for s in range(0, len(uniq), 20000):                        # chunk memory
        chunk = uniq[s:s + 20000]
        logp = base[None, :] - np.outer(chunk, g) / 50.0
        best[s:s + 20000] = logp.argmax(axis=1) + 1
    return best[inv]


# ----------------------------------------------------------------- data input
def parse_floats(s):
    return [float(x) for x in FLOAT_RE.findall(s)] if isinstance(s, str) else []


def load_segments(run_dir, it, check_align=True):
    prefix = f"{run_dir.rstrip('/')}/iter{it}"
    ibd = pd.read_csv(f"{prefix}.ibd.gz", sep=r"\s+", header=None,
                      names=["id1", "hap1", "id2", "hap2", "chr", "start", "end", "cM"])
    tm = pd.read_csv(f"{prefix}.tmrca.gz", sep="\t",
                     usecols=["chrom_index", "chromosome", "proportion", "tmrca"], dtype=str)
    if len(ibd) != len(tm):
        sys.exit(f"Row mismatch: {len(ibd)} IBD rows vs {len(tm)} TMRCA rows. "
                 "Rows are expected to correspond one-to-one.")

    if check_align:
        # tmrca.gz row i must describe ibd.gz row i.  concat_tmrca.py writes the chromosome
        # and the row index WITHIN that chromosome's IBD file, so both can be verified.
        try:
            tm_chr = tm["chromosome"].astype(int).values
            tm_idx = tm["chrom_index"].astype(int).values
        except ValueError:
            sys.exit(f"iter{it}.tmrca.gz has non-numeric chromosome/chrom_index values "
                     "(a stray header line?).")
        ibd_idx = ibd.groupby("chr").cumcount().values
        bad = np.flatnonzero((tm_chr != ibd["chr"].values) | (tm_idx != ibd_idx))
        if len(bad):
            i = bad[0]
            sys.exit(f"iter{it}: IBD and TMRCA rows are NOT aligned ({len(bad):,} rows differ; first "
                     f"at row {i}: ibd chr {ibd['chr'].iloc[i]} / position-in-chr {ibd_idx[i]} vs "
                     f"tmrca chromosome {tm_chr[i]} / chrom_index {tm_idx[i]}). "
                     "Pass --skip-align-check only if you know why this is expected.")

    props = tm["proportion"].map(parse_floats)
    times = tm["tmrca"].map(parse_floats)

    seg_t = np.full(len(tm), np.nan)
    for i, (p, t) in enumerate(zip(props, times)):
        if p and len(p) == len(t):
            seg_t[i] = t[int(np.argmax(p))]          # dominant MRCA node
    ibd["seg_tmrca"] = seg_t
    return ibd


def add_pair_tmrca(ibd):
    a = np.where(ibd["id1"] < ibd["id2"], ibd["id1"], ibd["id2"])
    b = np.where(ibd["id1"] < ibd["id2"], ibd["id2"], ibd["id1"])
    ibd["pair"] = pd.Series(a, index=ibd.index) + "|" + pd.Series(b, index=ibd.index)
    ibd["pair_tmrca"] = ibd.groupby("pair")["seg_tmrca"].transform("min")
    return ibd


# ------------------------------------------------------------------- plotting
# Panel A (reviewer-suggested view): segments binned by length, coloured by the
# segment's own TMRCA bin.  Bins follow the mock figure; the length bins are
# configurable (--len-bins) because the run's IBD files start at ~2 cM.
TMRCA_EDGES = [0, 1, 2, 5, 10, 20, 50, np.inf]
TMRCA_LABELS = ["1", "2", "3–5", "6–10", "11–20", "21–50", ">50"]
DEFAULT_LEN_EDGES = "2,3,5,10,20,inf"


def len_labels(edges):
    return [f"≥{lo:g}" if np.isinf(hi) else f"{lo:g}–{hi:g}"
            for lo, hi in zip(edges[:-1], edges[1:])]


def analytical_panel_a(N_long, len_edges):
    """Panel A from the model alone: P(segment TMRCA bin | segment length in [lo, hi)).

    Uses only the N(g) trajectory -- no simulated segments.  Eq. 1 is a conditional,
    P(g | l) = f(g, l) / gamma_0(l), of the joint
        f(g, l) = c(g) (g/50)^2 exp(-l g / 50),   c(g) = prod_{g'<g}(1 - 1/(2N[g'])) / (2N[g]),
    where gamma_0(l) = sum_g f(g, l) is the model's own (unnormalised) length density.
    For a length BIN the probability of TMRCA g is therefore the joint integrated over the
    bin, normalised over g:
        F(g; lo, hi) = int_lo^hi f(g, l) dl = c(g) (g/50) [exp(-lo g/50) - exp(-hi g/50)]
    i.e. lengths inside the bin are weighted by the model's own length distribution
    (short segments dominate), not uniformly.  Normalised over g = 1..len(N_long).
    Rows: the same length bins as the simulated panel A; columns: TMRCA bins.
    """
    g, base = _eq1_terms(N_long)
    log_g50 = np.log(g / 50.0)
    rows = {}
    for lab, lo, hi in zip(len_labels(len_edges), len_edges[:-1], len_edges[1:]):
        # log F = log c(g) + log(g/50) - lo g/50 + log(1 - exp(-(hi-lo) g/50));  hi = inf -> last term 0
        logF = (base - log_g50) - lo * g / 50.0 + np.log(-np.expm1(-(hi - lo) * g / 50.0))
        p = np.exp(logF - logF.max())
        p /= p.sum()
        rows[lab] = [p[(g > a) & (g <= b)].sum()
                     for a, b in zip(TMRCA_EDGES[:-1], TMRCA_EDGES[1:])]
    return pd.DataFrame(rows, index=TMRCA_LABELS).T


def make_groups(grouped, overflow=False):
    if grouped:
        groups = [("1", [1]), ("2", [2]), ("3", [3]), ("4–5", [4, 5]), ("6–10", list(range(6, 11)))]
    else:
        groups = [(str(t), [t]) for t in range(1, 11)]
    if overflow:                                    # pair TMRCA > 10, drawn in grey
        groups.append((">10", [SHOW_G + 1]))
    return groups


def legend_below(ax, title, rows="auto"):
    """Legend under the axes in 1 or 2 rows (rows="auto": 1 row if it fits the axes width).

    Entries read left to right in stacking order (first category = bottom of the stack).
    Matplotlib fills legend columns top-to-bottom, so a 2-row legend is re-ordered to keep
    row-major reading order (the last cell is padded with a blank entry when n is odd).
    """
    handles, labels = ax.get_legend_handles_labels()
    n = len(handles)
    title = title.replace("\n", " ")

    def build(nrows):
        ncol = math.ceil(n / nrows)
        hs, ls = list(handles), list(labels)
        if nrows == 2:
            pad = 2 * ncol - n
            hs += [Line2D([], [], color="none")] * pad
            ls += [""] * pad
            order = [i for j in range(ncol) for i in (j, j + ncol)]
            hs, ls = [hs[i] for i in order], [ls[i] for i in order]
        return ax.legend(hs, ls, title=title, ncol=ncol, frameon=False, loc="upper center",
                         bbox_to_anchor=(0.5, -0.17), borderaxespad=0, fontsize=9,
                         title_fontsize=9, columnspacing=1.2, handlelength=1.3,
                         handletextpad=0.5)

    if str(rows) in ("1", "2"):
        return build(int(rows))
    leg = build(1)
    ax.figure.canvas.draw()
    if leg.get_window_extent().width > ax.get_window_extent().width:
        leg.remove()
        leg = build(2)
    return leg


def draw_stacked(ax, props, groups, cmap_name, xlabel, legend_title,
                 xlabels=None, label_min=0.06, legend_rows="auto"):
    """Stacked bars. props: DataFrame, rows = bars, columns = categories, rows sum to 1."""
    cmap = plt.get_cmap(cmap_name)
    n_col = sum(1 for lab, _ in groups if lab != ">10")
    ramp = [cmap(x) for x in np.linspace(0.05, 0.95, n_col)]
    colors, k = [], 0
    for lab, _ in groups:                           # the ">10" overflow class is neutral grey
        if lab == ">10":
            colors.append((0.74, 0.74, 0.74, 1.0))
        else:
            colors.append(ramp[k])
            k += 1
    xs = np.arange(len(props)) if xlabels is not None else np.asarray(props.index)
    bottom = np.zeros(len(props))
    for (label, members), col in zip(groups, colors):
        vals = props[list(members)].sum(axis=1).values
        ax.bar(xs, vals, bottom=bottom, width=0.62, color=col,
               edgecolor="white", linewidth=0.6, label=label)
        for x, v, b in zip(xs, vals, bottom):
            if v >= label_min:
                dark = np.mean(col[:3]) < 0.5
                ax.text(x, b + v / 2, f"{v * 100:.0f}%", ha="center", va="center",
                        fontsize=7, color="white" if dark else "black")
        bottom += vals
    ax.set_xticks(xs)
    if xlabels is not None:
        ax.set_xticklabels(xlabels)
    ax.set_xlim(xs.min() - 0.6, xs.max() + 0.6)
    ax.set_ylim(0, 1)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_yticklabels([f"{int(t * 100)}%" for t in np.linspace(0, 1, 6)])
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Proportion of segments")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    legend_below(ax, legend_title, legend_rows)


def plot(props, groups, path, label_min=0.06, legend_rows="auto"):
    fig, ax = plt.subplots(figsize=(7.6, 5.2), dpi=300)
    draw_stacked(ax, props, groups, "viridis",
                 xlabel="Most likely segment age (generations)",
                 legend_title="Pair TMRCA (generations)", label_min=label_min,
                 legend_rows=legend_rows)
    fig.savefig(path, bbox_inches="tight")          # tight bbox also keeps the legend in frame
    plt.close(fig)


def plot_two_panel(a_props, b_props, b_groups, path, b_label_min=0.06, legend_rows="auto"):
    """A: segment TMRCA by segment length (reviewer's design).  B: Fig 2 breakdown."""
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(13.4, 5.4), dpi=300,
                                   gridspec_kw=dict(width_ratios=[1, 1.5], wspace=0.3))
    draw_stacked(axa, a_props, [(c, [c]) for c in a_props.columns], "plasma",
                 xlabel="IBD segment length (cM)", legend_title="Segment TMRCA (generations)",
                 xlabels=list(a_props.index), label_min=0.05, legend_rows=legend_rows)
    draw_stacked(axb, b_props, b_groups, "viridis",
                 xlabel="Most likely segment age (generations)",
                 legend_title="Pair TMRCA (generations)", label_min=b_label_min,
                 legend_rows=legend_rows)
    axa.text(-0.13, 1.05, "A", transform=axa.transAxes, fontsize=16, fontweight="bold")
    axb.text(-0.08, 1.05, "B", transform=axb.transAxes, fontsize=16, fontweight="bold")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------- main
def find_iters(run_dir):
    found = []
    for f in os.listdir(run_dir):
        m = re.fullmatch(r"iter(\d+)\.ibd\.gz", f)
        if not m:
            continue
        it = int(m.group(1))
        if os.path.exists(f"{run_dir.rstrip('/')}/iter{it}.tmrca.gz"):
            found.append(it)
        else:
            print(f"  skipping iter{it}: no iter{it}.tmrca.gz")
    return sorted(found)


def panel_a_counts(ibd, len_edges):
    """Rows: segment-length bins (cM).  Columns: segment-TMRCA bins (generations)."""
    d = ibd[ibd["seg_tmrca"].notna()]
    t = pd.cut(np.rint(d["seg_tmrca"]).clip(lower=1), TMRCA_EDGES,
               labels=TMRCA_LABELS, right=True)
    l = pd.cut(d["cM"], len_edges, labels=len_labels(len_edges), right=False)
    ct = pd.crosstab(l, t)
    return ct.reindex(index=len_labels(len_edges), columns=TMRCA_LABELS, fill_value=0).astype(int)


def panel_a_missing(ibd, len_edges):
    """Per length bin: all detected segments, and those with no TMRCA entry (dropped from A)."""
    labels = len_labels(len_edges)
    l = pd.cut(ibd["cM"], len_edges, labels=labels, right=False)
    d = ibd.assign(_len=l).dropna(subset=["_len"])
    g = d.groupby("_len", observed=False)
    out = pd.DataFrame({"segments": g.size(), "no_tmrca": g["seg_tmrca"].apply(lambda s: int(s.isna().sum()))})
    return out.reindex(labels).fillna(0).astype(int)


def counts_for_iter(run_dir, it, N, over10, len_edges, check_align=True):
    """Fig 2 table (age x pair TMRCA), simulation-based panel A tables, and stats.

    over10 controls segments whose pair TMRCA is > 10:
      "separate"  extra class SHOW_G+1 (">10"); bars cover ALL segments with a resolved pair
      "drop"      excluded, so bars are renormalised over pair TMRCA 1..10
      "clip"      folded into the 10 class (10 then means "10 or more")
    """
    ibd = add_pair_tmrca(load_segments(run_dir, it, check_align))
    ibd["age"] = most_likely_age(ibd["cM"].values, N)

    d = ibd[(ibd["age"] <= SHOW_G) & ibd["pair_tmrca"].notna()].copy()
    d["pt"] = np.rint(d["pair_tmrca"]).astype(int).clip(lower=1)
    n_over = int((d["pt"] > SHOW_G).sum())
    cols = list(range(1, SHOW_G + 1))
    if over10 == "drop":
        d = d[d["pt"] <= SHOW_G]
    elif over10 == "clip":
        d = d.assign(pt=d["pt"].clip(upper=SHOW_G))
    else:
        d = d.assign(pt=d["pt"].clip(upper=SHOW_G + 1))
        cols.append(SHOW_G + 1)

    counts = (d.groupby(["age", "pt"]).size().unstack(fill_value=0)
              .reindex(index=range(1, SHOW_G + 1), columns=cols, fill_value=0))
    a_counts = panel_a_counts(ibd, len_edges)
    a_missing = panel_a_missing(ibd, len_edges)
    stats = dict(segments=len(ibd), pairs=ibd["pair"].nunique(),
                 no_tmrca=int(ibd["seg_tmrca"].isna().sum()),
                 unresolved_pair=int(ibd["pair_tmrca"].isna().sum()),
                 age_ok=int(((ibd["age"] <= SHOW_G) & ibd["pair_tmrca"].notna()).sum()),
                 over=n_over, a_used=int(a_counts.values.sum()))
    return counts, a_counts, a_missing, stats


def to_props(counts):
    return counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="run directory with iter{N}.* files")
    ap.add_argument("--iter", type=int, nargs="*", default=None,
                    help="iteration(s) to use (default: all found in the run directory)")
    ap.add_argument("--demo-file", required=True, help="python file defining the demography")
    ap.add_argument("--demo-object", default="ooa2")
    ap.add_argument("--pop", default="pop_0", help="sampled population (IBD-sims samples pop_0)")
    ap.add_argument("--over10", choices=["separate", "drop", "clip"], default="separate",
                    help="segments whose pair TMRCA is > 10: show as a separate grey '>10' class "
                         "(default; bars cover ALL segments), drop them (bars renormalised over "
                         "pair TMRCA 1-10), or clip them into the 10 class")
    ap.add_argument("--two-panel", action="store_true",
                    help="also write a two-panel figure: A = segment TMRCA by segment length "
                         "(reviewer's design), B = the Fig 2 breakdown")
    ap.add_argument("--a-analytical", action="store_true",
                    help="plot panel A from Eq. 1 and the N(g) trajectory alone (integrated over each "
                         "--len-bins bin) instead of from the simulated segments' true TMRCAs")
    ap.add_argument("--a-max-g", type=int, default=1000,
                    help="analytical panel A: normalise Eq. 1 over g = 1..this (default 1000; the "
                         "paper's age assignment uses 1..100, but for short segments a sizeable "
                         "share of the distribution lies beyond 100 generations)")
    ap.add_argument("--panel-b", choices=["full", "grouped"], default="full",
                    help="which Fig 2 version to use as panel B (default: full, 10 classes)")
    ap.add_argument("--len-bins", default=DEFAULT_LEN_EDGES,
                    help="comma-separated segment-length bin edges in cM for panel A "
                         f"(default: {DEFAULT_LEN_EDGES}); bins are [lo, hi)")
    ap.add_argument("--legend-rows", choices=["auto", "1", "2"], default="auto",
                    help="legends sit below each x-axis in 1 or 2 rows (default auto: 1 row if it "
                         "fits the panel width, else 2)")
    ap.add_argument("--skip-align-check", action="store_true",
                    help="do not verify that tmrca.gz rows line up with ibd.gz rows "
                         "(chromosome and within-chromosome index)")
    ap.add_argument("--out", default="fig2_new", help="output prefix")
    args = ap.parse_args()

    len_edges = [float(x) for x in args.len_bins.split(",")]

    N_long = load_ne_trajectory(args.demo_file, args.demo_object, args.pop,
                                max_g=max(MAX_G, args.a_max_g))
    N = N_long[:MAX_G]                      # g = 1..100: the paper's age-assignment range
    N_a = N_long[:args.a_max_g]             # g = 1..a_max_g: range for theory panel A
    print(f"N(g) for {args.pop}: g=1 -> {N[0]:.0f}, g=10 -> {N[9]:.0f}, g=100 -> {N[99]:.0f}"
          + (f", g={len(N_a)} -> {N_a[-1]:.0f}" if len(N_a) > MAX_G else ""))

    iters = args.iter if args.iter else find_iters(args.run)
    if not iters:
        sys.exit(f"No iter*.ibd.gz with a matching tmrca.gz found in {args.run}")
    print(f"iterations used ({len(iters)}): {iters if len(iters) <= 12 else str(iters[:12])[:-1] + ', ...]'}")

    per_iter, per_iter_a, per_iter_miss, tot = {}, {}, {}, {}
    for it in iters:
        counts, a_counts, a_miss, st = counts_for_iter(
            args.run, it, N, args.over10, len_edges, not args.skip_align_check)
        per_iter[it], per_iter_a[it], per_iter_miss[it] = counts, a_counts, a_miss
        for k, v in st.items():
            tot[k] = tot.get(k, 0) + v
        print(f"  iter{it}: {st['segments']:,} segments, {st['pairs']:,} pairs")

    counts = sum(per_iter.values())
    props = to_props(counts)
    a_counts = sum(per_iter_a.values())
    a_props = to_props(a_counts)
    a_miss = sum(per_iter_miss.values())
    a_miss["share_without_tmrca"] = a_miss["no_tmrca"] / a_miss["segments"].replace(0, np.nan)
    analytical = analytical_panel_a(N_a, len_edges)
    a_plot = analytical if args.a_analytical else a_props

    print(f"\ntotal segments: {tot['segments']:,}   (pairs summed over iterations: {tot['pairs']:,})")
    print(f"segments with no TMRCA row entry: {tot['no_tmrca'] / tot['segments']:.1%}")
    print(f"segments whose pair has no resolved TMRCA: {tot['unresolved_pair'] / tot['segments']:.1%}")
    print(f"segments with most likely age <= {SHOW_G} and resolved pair: {tot['age_ok']:,}")
    how = {"separate": "shown as a separate '>10' class", "drop": "DROPPED: bars are renormalised "
           "over pair TMRCA 1-10", "clip": "clipped into the 10 class"}[args.over10]
    print(f"  of these, pair TMRCA > {SHOW_G}: {tot['over'] / max(tot['age_ok'], 1):.1%}  ({how})")
    print(f"segments used in simulation-based panel A (have a TMRCA, length within bins): {tot['a_used']:,}")

    # outputs: pooled counts/proportions, plus per-iteration long table
    counts.to_csv(f"{args.out}_counts.csv")
    props.to_csv(f"{args.out}_proportions.csv")
    a_counts.to_csv(f"{args.out}_A_counts.csv")
    a_props.to_csv(f"{args.out}_A_proportions.csv")
    a_miss.to_csv(f"{args.out}_A_missing_tmrca.csv")
    analytical.to_csv(f"{args.out}_A_analytical_proportions.csv")
    long = pd.concat({it: to_props(c).stack() for it, c in per_iter.items()},
                     names=["iter", "age", "pair_tmrca"])
    long.rename("proportion").reset_index().to_csv(f"{args.out}_proportions_by_iter.csv", index=False)
    sd = None
    if len(iters) > 1:
        sd = long.groupby(level=["age", "pair_tmrca"]).std().unstack()
        sd.to_csv(f"{args.out}_proportions_sd_across_iters.csv")

    pd.set_option("display.float_format", lambda x: f"{x * 100:5.1f}")
    print("\nPooled proportion of segments (%), rows = most likely age, cols = pair TMRCA:")
    print(props.to_string())
    cum3 = props[[1, 2, 3]].sum(axis=1)
    print("\nCumulative share from pair TMRCA <= 3 (~1st-4th degree), by age:")
    print(cum3.to_string())

    print("\nCheck against the published Fig 2 (approximate values read off the figure):")
    pub = {(1, 1): 65, (1, 2): 16, (1, 3): 16, (2, 1): 36, (2, 2): 34, (2, 3): 10,
           (3, 2): 24, (6, 1): 6, (10, 10): 12}
    for (g, t), v in pub.items():
        s = f"  (SD across iters {sd.loc[g, t] * 100:4.1f})" if sd is not None else ""
        print(f"  age {g:>2}, pair TMRCA {t:>2}:  published ~{v:>2}%   here {props.loc[g, t] * 100:5.1f}%{s}")
    print(f"  age  1, TMRCA<=3: published ~97%   here {cum3.loc[1] * 100:5.1f}%")
    print(f"  age  8, TMRCA<=3: published ~17%   here {cum3.loc[8] * 100:5.1f}%")

    print("\nPanel A, SIMULATION (true segment TMRCAs from tmrca.gz), %: rows = length bin (cM)"
          + ("" if args.a_analytical else "   [plotted]"))
    print(a_props.to_string())
    print("(counts per length bin: " + ", ".join(f"{k}: {int(v):,}" for k, v in a_counts.sum(axis=1).items()) + ")")
    print("\nSimulation-based panel A: segments dropped for having no TMRCA entry, by length bin "
          "(share = %):")
    print(a_miss.to_string())
    print("(if this share is large or uneven across bins, the simulation-based panel A describes only "
          "the subset of segments whose true MRCA could be resolved)")

    print(f"\nPanel A, ANALYTICAL (Eq. 1 integrated over each length bin; N(g) only; normalised over "
          f"g=1..{len(N_a)}), %" + ("   [plotted]" if args.a_analytical else ""))
    print(analytical.to_string())
    print("\nSimulation minus analytical, percentage points (rows = length bin):")
    print((a_props - analytical).to_string())

    overflow = args.over10 == "separate"
    plot(props, make_groups(False, overflow), f"{args.out}_full.png", legend_rows=args.legend_rows)
    plot(props, make_groups(True, overflow), f"{args.out}_grouped.png", label_min=0.04,
         legend_rows=args.legend_rows)
    written = [f"{args.out}_full.png", f"{args.out}_grouped.png"]
    if args.two_panel:
        grouped = args.panel_b == "grouped"
        plot_two_panel(a_plot, props, make_groups(grouped, overflow), f"{args.out}_twopanel.png",
                       b_label_min=0.04 if grouped else 0.06, legend_rows=args.legend_rows)
        written.append(f"{args.out}_twopanel.png")
    print(f"\nWrote {', '.join(written)} and CSVs "
          f"(*_counts, *_proportions, *_A_analytical_proportions, *_A_counts, *_A_proportions, "
          f"*_A_missing_tmrca, *_proportions_by_iter"
          + (", *_proportions_sd_across_iters)" if sd is not None else ")"))


if __name__ == "__main__":
    main()
