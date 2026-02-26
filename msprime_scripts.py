import msprime
import pandas as pd
import numpy as np
import itertools as it
import argparse

from write_vcf import write_vcf

def parse_arguments():
    # Create the parser
    parser = argparse.ArgumentParser(description='msprime simulation script')

    # Which model to use
    parser.add_argument("--mode",
                        choices=["coalescent", "pedigree"],
                        help="Simulation type")

    parser.add_argument("--dtwf_end", type=int, default=25)

    parser.add_argument("--constant_Ne", type=int, default=10_000)

    parser.add_argument("--samples", type=int, default=1000)

    parser.add_argument("--pedigree", type=str)

    parser.add_argument("--output", type=str)

    parser.add_argument("--chr", type=int, default=1)
    
    return parser.parse_args()

def simulate_coalescent(samples, Ne, recomb_rate=1e-8, seq_lenth=100_000_000, end_time=None):
    ts = msprime.sim_ancestry(
        samples=samples,
        sequence_length=seq_lenth,
        recombination_rate=recomb_rate,
        population_size=Ne,
        discrete_genome=True,
        end_time=end_time,
        initial_state=None
    )

    return ts

def simulate_pedigree(pedigree_file, end_time, Ne, sequence_length=100_000_000, recomb_rate=1e-8):

    pedigree = msprime.parse_pedigree(open(pedigree_file), sequence_length=sequence_length)

    ped_ts = msprime.sim_ancestry(initial_state=pedigree,
                                model="fixed_pedigree",
                                recombination_rate=recomb_rate,
                                end_time=end_time,
                                discrete_genome=True)

    ts = msprime.sim_ancestry(initial_state=ped_ts,
                             recombination_rate=recomb_rate,
                             population_size=Ne,
                             discrete_genome=True)

    return ts

def simulate_chromosome(args, mutation_rate=1e-8):

    if args.mode == "coalescent":

        ts = simulate_coalescent(args.samples,
                                 args.constant_Ne)

    elif args.mode == "pedigree":

        ts = simulate_pedigree(args.pedigree,
                               args.dtwf_end,
                               args.constant_Ne)

    name = f"Ne{args.constant_Ne}_n{args.samples}_{args.mode}"

    ts.dump(f"{args.output}_{name}_chr{args.chr}.tree")

    mut_ts = msprime.sim_mutations(ts, rate=mutation_rate)

    write_vcf(mut_ts, f"{args.output}_{name}_chr{args.chr}", args.chr)


if __name__ == "__main__":

    args = parse_arguments()

    simulate_chromosome(args)