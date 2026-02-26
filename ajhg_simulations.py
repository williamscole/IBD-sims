import multiprocessing as mp
import pandas as pd
import numpy as np
import msprime
from functools import partial
import itertools as it
from write_vcf import write_vcf
import sys

PARAMS = {
    'Ne': 10_000,
    'samples': 1_000,
    'recomb_rate': 1e-8,
    'sequence_length': 100_000_000,
    'end_time': 100
}

def simulate_chromosome(chrom, params):
    """Simulate a single chromosome and get its IBD segments"""
    ts = msprime.sim_ancestry(
        samples=params['samples'],
        sequence_length=params['sequence_length'],
        recombination_rate=params['recomb_rate'],
        population_size=params['Ne'],
        discrete_genome=True,
    )

    ts.dump(f"trees/ajhg_sim_chr{chrom}.tree")

    mut_ts = msprime.sim_mutations(ts, rate=10e-8)

    write_vcf(mut_ts, f"vcfs/ajhg_sim_chr{chrom}")


if __name__ == "__main__":

    chrom = sys.argv[-1]

    simulate_chromosome(chrom, PARAMS)