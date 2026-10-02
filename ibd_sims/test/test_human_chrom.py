#!/usr/bin/env python
"""Diagnose a human-recombination-map (end_chr: 22) run, one stage at a time.

Lives in <repo>/ibd_sims/test/ and can be run from anywhere. Each stage calls the repo's own code and reports time and
peak memory, and the script stops at the first failure so you can see where the
real pipeline breaks.

    python ibd_sims/test/test_human_chrom.py OUT_DIR --chrom 22 --samples 100
    python ibd_sims/test/test_human_chrom.py OUT_DIR --chrom 1  --samples 1000   # scale test
    python ibd_sims/test/test_human_chrom.py OUT_DIR --check-all                 # only check the 22 map files
    python ibd_sims/test/test_human_chrom.py OUT_DIR --chrom 22 --samples 100 --hapibd        # also VCF + hap-ibd, compared with tskit
    python ibd_sims/test/test_human_chrom.py OUT_DIR --chrom 22 --samples 100 --hapibd-only   # VCF + hap-ibd only

The chr1 map path is read from hapmap_chr1 in <repo>/setup.yaml (override with --hapmap);
other chromosomes are found by replacing "chr1" with "chr{n}", exactly as the pipeline does.

Stages: 1 map files | 2 genome setup | 3 simulate | 4 tskit IBD + map | 5 TMRCA |
        6/7 (optional, --hapibd) VCF writing + hap-ibd, compared with the tskit calls
        --hapibd-only skips 4-5 and the comparison, and adds 8 (TMRCA on the hap-ibd calls)
"""
import argparse
import os
import resource
import sys
import time
import traceback
from pathlib import Path

# this file lives in <repo>/ibd_sims/test/, so the repo root is two levels up
REPO = str(Path(__file__).resolve().parents[2])


def mem_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6   # Linux: KB -> GB


class Stage:
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.t = time.time()
        print(f"\n=== {self.name}", flush=True)
        return self

    def __exit__(self, etype, exc, tb):
        dt = time.time() - self.t
        if etype is None:
            print(f"--- PASS  ({dt:.1f}s, peak memory so far {mem_gb():.2f} GB)", flush=True)
            return False
        print(f"--- FAIL  ({dt:.1f}s): {etype.__name__}: {exc}", flush=True)
        traceback.print_exception(etype, exc, tb, limit=-4)
        sys.exit(1)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("out_dir")
    p.add_argument("--hapmap", default=None, help="chr1 HapMap path (default: hapmap_chr1 from <repo>/setup.yaml)")
    p.add_argument("--repo", default=REPO, help="repo root (default: inferred from this file's location)")
    p.add_argument("--chrom", type=int, default=22)
    p.add_argument("--samples", type=int, default=100, help="diploid individuals")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--demo-path", default=None, help="demography file (default: <repo>/ibd_sims/demography.py)")
    p.add_argument("--demo-object", default="constant_Ne", help="demography object name (default constant_Ne)")
    p.add_argument("--check-all", action="store_true", help="only check that all 22 map files exist and parse")
    p.add_argument("--hapibd", action="store_true",
                   help="also do the VCF + hap-ibd route, using hap_ibd_jar and maf_pickle from setup.yaml, "
                        "and compare segment counts with the tskit route")
    p.add_argument("--hapibd-only", action="store_true",
                   help="the VCF + hap-ibd route only: skip the tskit IBD stages and the comparison "
                        "(what the pipeline does with tskit_ibd: false)")
    p.add_argument("--vcf", action="store_true", help="only the VCF stage (no hap-ibd)")
    p.add_argument("--snps-pkl", default=None, help="SNP pickle (default: maf_pickle from setup.yaml)")
    p.add_argument("--hap-ibd-jar", default=None, help="hap-ibd jar (default: hap_ibd_jar from setup.yaml)")
    p.add_argument("--gb", type=int, default=8, help="memory (GB) to give hap-ibd")
    a = p.parse_args()

    sys.path.insert(0, os.path.join(a.repo, "ibd_sims"))
    os.makedirs(a.out_dir, exist_ok=True)
    chrom = a.chrom
    import pandas as pd
    import yaml
    setup = yaml.safe_load(open(os.path.join(a.repo, "setup.yaml")))
    if a.hapmap is None:
        a.hapmap = setup["hapmap_chr1"]
    if a.hapibd or a.hapibd_only:
        a.vcf = True
        a.hap_ibd_jar = a.hap_ibd_jar or setup["hap_ibd_jar"]
        a.snps_pkl = a.snps_pkl or setup["maf_pickle"]
    print(f"HapMap (chr1 path): {a.hapmap}")

    # 1 ── the map files ------------------------------------------------------------
    chroms = range(1, 23) if a.check_all else [chrom]
    with Stage(f"1. HapMap files ({'all 22' if a.check_all else f'chr{chrom}'})"):
        for c in chroms:
            f = a.hapmap.replace("chr1", f"chr{c}")
            if not os.path.exists(f):
                raise FileNotFoundError(f"{f}  (made by replacing 'chr1' with 'chr{c}' in the whole path)")
            df = pd.read_csv(f, sep=r"\s+")
            need = {"Position(bp)", "Rate(cM/Mb)", "Map(cM)"}
            if not need <= set(df.columns):
                raise ValueError(f"{f}: columns {list(df.columns)}, expected {sorted(need)}")
            print(f"  chr{c}: {len(df):,} rows, last position {int(df['Position(bp)'].iloc[-1]):,}, "
                  f"genetic length {df['Map(cM)'].iloc[-1] - df['Map(cM)'].iloc[0]:.1f} cM")
    if a.check_all:
        print("\nAll map files look fine.")
        return

    # 2 ── genome setup (what GenomeSetup does for end_chr: 22) -----------------------
    import numpy as np
    import msprime
    from simulations import GenomeSetup, Simulation
    args = {"end_chr": 22, "hapmap_chr1": a.hapmap, "samples": a.samples, "seed": a.seed,
            "iter_n": 1, "chrom": chrom,
            "custom_sim": {"path": None, "object": None},
            "custom_demo": {"path": a.demo_path or os.path.join(a.repo, "ibd_sims", "demography.py"),
                            "object": a.demo_object},
            "pedigree": {"pedigree_mode": False}}
    with Stage("2. GenomeSetup.create (read_hapmap -> RateMap)"):
        seq_len, rate = GenomeSetup.create(args, chrom)
        print(f"  sequence length {seq_len:,} bp; rate type {type(rate).__name__}")
        if hasattr(rate, "rate"):
            print(f"  intervals {len(rate.rate):,}, missing (NaN) intervals {int(np.isnan(rate.rate).sum())}, "
                  f"total {rate.get_cumulative_mass(seq_len) * 100:.1f} cM")

    # 3 ── simulate with the repo's own Simulation.create (coalescent, constant Ne) ----
    with Stage(f"3. Simulation.create: {a.samples} diploids on chr{chrom}"):
        ts, rate, demog, seed = Simulation.create(args, chrom, os.path.join(a.out_dir, "iter1"))
        print(f"  {ts.num_samples} haplotypes, {ts.num_trees:,} trees, {ts.num_sites:,} sites, "
              f"{ts.num_edges:,} edges")

    # 4 ── tskit IBD + map (what sim() does when tskit_ibd: true) ----------------------
    from simulations import add_tmrca
    prefix = os.path.join(a.out_dir, f"iter1_chr{chrom}")
    df = None
    if not a.hapibd_only:
        from ibd_from_ts import ts_ibd_pipeline
        with Stage("4. ts_ibd_pipeline (tskit IBD, .map, .samples)"):
            df, _ = ts_ibd_pipeline(ts, a.out_dir, 1, chrom, rate)
            m = pd.read_csv(f"{prefix}.map", sep=r"\s+", header=None, names=["chrom", "id", "cm", "bp"])
            print(f"  {len(df):,} IBD segments, mean {df.cM.mean() if len(df) else 0:.2f} cM")
            print(f"  map: {len(m):,} rows, last bp {m.bp.iloc[-1]:,} (sequence length {seq_len:,}), "
                  f"total {m.cm.iloc[-1]:.1f} cM")
            assert m.bp.iloc[-1] == seq_len, "map does not end at the sequence length"

        # 5 ── TMRCA (reads the .ibd.gz back) -----------------------------------------
        with Stage("5. add_tmrca"):
            add_tmrca(prefix, ts, False)
            print(f"  wrote {prefix}.tmrca.pkl")

    # 6/7 ── hap-ibd route -------------------------------------------------------------
    if a.vcf or a.hap_ibd_jar:
        from write_vcf import write_vcf
        pkl = a.snps_pkl or setup["maf_pickle"]
        if not os.path.isabs(pkl) and not os.path.exists(pkl):     # relative paths: try the repo root
            pkl = os.path.join(a.repo, pkl)
        hprefix = prefix if a.hapibd_only else prefix + "_hapibd"
        with Stage("6. write_vcf (SNP thinning, VCF, hap-ibd map)"):
            write_vcf(ts, hprefix, chrom, rate, a.seed, snps_pkl=pkl)
            print(f"  VCF {os.path.getsize(hprefix + '.vcf.gz') / 1e6:.1f} MB")
    if a.hap_ibd_jar:
        from simulations import run_hapibd
        with Stage(f"7. hap-ibd ({a.gb} GB)"):
            ok = run_hapibd(hprefix, a.gb, hapibd_jar=a.hap_ibd_jar)
            if not ok or not os.path.exists(hprefix + ".ibd.gz"):
                raise RuntimeError("hap-ibd failed (see its error output above)")
            n = sum(1 for _ in __import__("gzip").open(hprefix + ".ibd.gz"))
            print(f"  hap-ibd called {n:,} segments" + (f" (tskit route: {len(df):,})" if df is not None else ""))
        if a.hapibd_only:
            with Stage("8. add_tmrca on the hap-ibd calls"):
                add_tmrca(hprefix, ts, False)
                print(f"  wrote {hprefix}.tmrca.pkl")

    print(f"\nAll stages passed. Peak memory {mem_gb():.2f} GB.")


if __name__ == "__main__":
    main()
