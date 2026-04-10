import pandas as pd
import numpy as np
import pickle as pkl
import sys
import msprime
from simulations import GenomeSetup

CHROM22 = "/gpfs/data/sramacha/datasets/quebec_sims/simulated_chrom_22.ts.tsz"

PATH = "quebec_trees"

def simplify_quebec(chrom, args, seed=42):
    import tszip

    ts = tszip.decompress(CHROM22.replace("_22.", f"_{chrom}."))

    n_samples = np.where(ts.nodes_time==0)[0].shape[0] // 2

    np.random.seed(42)

    random_nodes = np.random.choice(np.arange(n_samples), args.get("samples", 1000), replace=False) * 2

    random_nodes = np.sort(np.concatenate((random_nodes, random_nodes+1)))

    ts = ts.simplify(samples=random_nodes)
    
    ts.dump(f"{PATH}/chr{chrom}_random_10000.trees")

def load_random_10000(chrom, args):

    iteration_seed = args.get("iteration_seed", 1000)

    np.random.seed(iteration_seed)

    ts = msprime.load(f"{PATH}/chr{chrom}_random_10000.trees")

    n_samples = np.where(ts.nodes_time==0)[0].shape[0] // 2

    random_nodes = np.random.choice(np.arange(n_samples), args.get("samples", 1000), replace=False) * 2
    random_nodes = np.sort(np.concatenate((random_nodes, random_nodes+1)))

    if len(random_nodes) < (n_samples * 2):
        ts = ts.simplify(samples=random_nodes)

    _, rate = GenomeSetup.create(args, chrom)

    return ts, rate


if __name__ == "__main__":

    for i in range(1, 23):
        simplify_quebec(i, {"samples": 10000})

