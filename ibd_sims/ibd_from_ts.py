import argparse
import gzip
import numbers
import os
import time

import msprime
import numpy as np
import pandas as pd


def test_ts(samples=50, seed=42):
    """Simple single-chromosome ts: constant Ne=10,000, 30 Mb, rate 1e-8.

    `samples` is the number of diploid individuals (haplotypes = 2 * samples).
    """
    ts = msprime.sim_ancestry(
        samples=samples,
        population_size=10_000,      # constant Ne
        sequence_length=30_000_000,  # 30 Mb "chromosome"
        recombination_rate=1e-8,
        discrete_genome=True,
        random_seed=seed,
    )

    return ts


def sample_node_map(ts):
    """Map sample node -> (j, hap), matching write_vcf's tsk_j convention.

    Individuals are ordered by their first sample node and named tsk_0..tsk_{N-1};
    hap is 1 or 2 by position within the individual's nodes.
    """
    sample_set = set(ts.samples().tolist())
    first = {}
    for n in ts.samples():
        ind = ts.node(n).individual
        if ind >= 0 and ind not in first:
            first[ind] = int(n)
    ordered = sorted(first, key=first.get)

    mapping = {}
    for j, ind in enumerate(ordered):
        nodes = [int(n) for n in ts.individual(ind).nodes if int(n) in sample_set]
        for h, n in enumerate(nodes, start=1):
            mapping[n] = (j, h)

    # Downstream code (get_node) assumes node = 2j + (hap == 2). Fail loudly if not.
    for n, (j, h) in mapping.items():
        if n != 2 * j + (h - 1):
            raise ValueError(
                f"Sample node {n} maps to (tsk_{j}, hap {h}), which breaks the "
                "node = 2j + (hap == 2) assumption used by get_node()."
            )
    return mapping


def make_cm_fn(rate=1e-8):
    """Return f(bp) -> cM.

    `rate` may be:
      - a number: constant recombination rate per bp per generation (default 1e-8)
      - an msprime.RateMap
      - the (sequence_length, rate) tuple returned by GenomeSetup.create(args, chrom)
    """
    if isinstance(rate, tuple):          # GenomeSetup.create output
        rate = rate[1]
    if isinstance(rate, numbers.Real):
        return lambda bp: np.asarray(bp, dtype=float) * float(rate) * 100
    if hasattr(rate, "get_cumulative_mass"):
        return lambda bp: rate.get_cumulative_mass(np.asarray(bp, dtype=float)) * 100
    raise TypeError(f"Unsupported rate type: {type(rate).__name__}")


def merge_segments(segs, cm_fn, gap_cm=0.0):
    """Merge a pair's segments that touch (or are within gap_cm of each other)."""
    segs = sorted(segs, key=lambda s: s.left)
    merged = []
    for s in segs:
        if merged:
            m = merged[-1]
            touching = s.left <= m["right"]
            near = gap_cm > 0 and (cm_fn(s.left) - cm_fn(m["right"])) <= gap_cm
            if touching or near:
                m["right"] = max(m["right"], s.right)
                m["pieces"].append((s.left, s.right, s.node))
                continue
        merged.append({"left": s.left, "right": s.right,
                       "pieces": [(s.left, s.right, s.node)]})
    return merged


def ibd_from_ts(ts, chrom, rate=1e-8, min_cm=2.0, gap_cm=0.0, min_span=0, max_time=None):
    """Call IBD directly from a tree sequence, in hap-ibd output layout.

    `rate` defaults to a constant 1e-8; it also accepts an msprime.RateMap or the
    (sequence_length, rate) tuple from GenomeSetup.create (see make_cm_fn).

    Returns a DataFrame with columns
        id1 hap1 id2 hap2 chr start end cM pieces
    where `pieces` is the list of (left, right, mrca_node) behind each merged
    segment (kept for TMRCA; not written to the .ibd.gz).
    """
    node_map = sample_node_map(ts)
    cm_fn = make_cm_fn(rate)

    kwargs = dict(min_span=min_span, store_segments=True)
    if max_time is not None:
        kwargs["max_time"] = max_time
    segments = ts.ibd_segments(**kwargs)

    rows = []
    for (n1, n2), seg_list in segments.items():
        j1, h1 = node_map[n1]
        j2, h2 = node_map[n2]
        if j1 == j2:          # same individual -> HBD, goes to .hbd.gz in hap-ibd
            continue
        for m in merge_segments(seg_list, cm_fn, gap_cm):
            length_cm = float(cm_fn(m["right"]) - cm_fn(m["left"]))
            if length_cm < min_cm:
                continue
            rows.append([f"tsk_{j1}", h1, f"tsk_{j2}", h2, str(chrom),
                         int(m["left"]), int(m["right"]), length_cm, m["pieces"]])

    return pd.DataFrame(rows, columns=["id1", "hap1", "id2", "hap2", "chr",
                                       "start", "end", "cM", "pieces"])


def write_ibd(df, prefix):
    """Write {prefix}.ibd.gz: 8 whitespace-separated columns, no header."""
    with gzip.open(f"{prefix}.ibd.gz", "wt") as f:
        for r in df.itertuples(index=False):
            f.write(f"{r.id1}\t{r.hap1}\t{r.id2}\t{r.hap2}\t{r.chr}\t"
                    f"{r.start}\t{r.end}\t{r.cM:.3f}\n")


def write_empty_hbd(prefix):
    """Write an empty {prefix}.hbd.gz so the pipeline sees the same files as hap-ibd.

    hap-ibd puts within-individual (HBD) segments here. Nothing downstream uses
    them, so this mode just creates the empty file.
    """
    with gzip.open(f"{prefix}.hbd.gz", "wt"):
        pass


def write_samples(ts, path):
    """Write the iter-level .samples file (tsk_0..tsk_{N-1}), as write_vcf does."""
    n = len({j for j, _ in sample_node_map(ts).values()})
    with open(path, "w") as f:
        for j in range(n):
            f.write(f"tsk_{j}\n")


def ts_ibd_pipeline(ts, out_dir, iter_n=1, chrom=1, rate=1e-8, **ibd_kwargs):
    """ts -> IBD calls -> {out_dir}/iter{n}_chr{chrom}.ibd.gz and iter{n}.samples.

    File names follow the repo's convention. Extra kwargs (min_cm, gap_cm,
    min_span, max_time) go to ibd_from_ts. Returns the calls DataFrame.
    """
    os.makedirs(out_dir, exist_ok=True)
    prefix = os.path.join(out_dir, f"iter{iter_n}_chr{chrom}")

    df = ibd_from_ts(ts, chrom, rate=rate, **ibd_kwargs)
    write_ibd(df, prefix)
    write_empty_hbd(prefix)
    write_samples(ts, os.path.join(out_dir, f"iter{iter_n}.samples"))
    return df, prefix


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Simulate a test chromosome (constant Ne=10,000) and call IBD "
                    "directly from the tree sequence, writing hap-ibd-style output.")
    p.add_argument("out_dir", help="directory to write outputs to (created if missing)")
    p.add_argument("--samples", type=int, default=50, help="diploid individuals (default 50)")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--iter", type=int, default=1, dest="iter_n", help="iteration number for file names")
    p.add_argument("--chrom", type=int, default=1)
    p.add_argument("--rate", type=float, default=1e-8, help="recombination rate per bp per gen")
    p.add_argument("--min-cm", type=float, default=2.0, help="minimum segment length in cM")
    p.add_argument("--gap-cm", type=float, default=0.0, help="merge segments separated by <= this many cM")
    p.add_argument("--min-span", type=int, default=20_000,
                   help="tskit min_span in bp; performance knob, lower = slower and more memory")
    p.add_argument("--max-time", type=float, default=None,
                   help="ignore ancestry older than this many generations")
    a = p.parse_args(argv)

    t0 = time.time()
    ts = test_ts(samples=a.samples, seed=a.seed)
    print(f"Simulated {ts.num_samples} haplotypes, {ts.num_trees} trees "
          f"({time.time() - t0:.1f}s)")

    t0 = time.time()
    df, prefix = ts_ibd_pipeline(ts, a.out_dir, iter_n=a.iter_n, chrom=a.chrom, rate=a.rate,
                              min_cm=a.min_cm, gap_cm=a.gap_cm,
                              min_span=a.min_span, max_time=a.max_time)
    print(f"Called {len(df)} IBD segments >= {a.min_cm} cM ({time.time() - t0:.1f}s)")
    if len(df):
        print(f"  mean length {df.cM.mean():.2f} cM, max {df.cM.max():.2f} cM")
    print(f"Wrote {prefix}.ibd.gz and {prefix}.hbd.gz (empty)")
    print(f"Wrote {os.path.join(a.out_dir, f'iter{a.iter_n}.samples')}")


if __name__ == "__main__":
    main()