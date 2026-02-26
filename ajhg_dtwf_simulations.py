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
    'end_time': 30
}

def simulate_chromosome(chrom, params):

    pedigree = msprime.parse_pedigree(open("wf_n1000_Ne10000.fam"), sequence_length=PARAMS["sequence_length"])

    ped_ts = msprime.sim_ancestry(initial_state=pedigree,
                                model="fixed_pedigree",
                                recombination_rate=PARAMS["recomb_rate"],
                                end_time=PARAMS["end_time"],
                                discrete_genome=True)

    """Simulate a single chromosome and get its IBD segments"""
    ts = msprime.sim_ancestry(initial_state=ped_ts,
        recombination_rate=params['recomb_rate'],
        population_size=params['Ne'],
        discrete_genome=True,
    )

    ts.dump(f"trees/ajhg_dtwf_sim_chr{chrom}.tree")

    mut_ts = msprime.sim_mutations(ts, rate=10e-8)

    write_vcf(mut_ts, f"vcfs/ajhg_dtwf_sim_chr{chrom}")

if __name__ == "__main__":

    chrom = sys.argv[-1]

    simulate_chromosome(chrom, PARAMS)

