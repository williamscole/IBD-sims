import yaml
import sys
import msprime
import numpy as np
import subprocess
import os
import pandas as pd
import importlib.util
import random
import hashlib
import itertools as it
import pickle as pkl
from pathlib import Path

from write_vcf import write_vcf
from wf_pedigree import create_pedigree
from maf_buckets import SNPs, BUCKETS


import warnings

pd.options.mode.chained_assignment = None

def load_config():
    config_path = Path(__file__).parent / "setup.yaml"
    return yaml.safe_load(open(config_path))

class WFSetup:
    @staticmethod
    def create(args, path, demography):

        debug = demography.debug()

        samples = args["samples"]
        gen_end = args["pedigree"]["gen_end"]
        mating = args["pedigree"]["mating"]
        iter_n = args["iter_n"]
        seed = args["iteration_seed"]

        return create_pedigree(samples, debug, f"{path}/iter{iter_n}", gen_end, mating, seed)



class GenomeSetup:
    @staticmethod
    def create(args, chrom):
        if args["end_chr"] == 30:
            sequence_length = 100_000_000
            rate = 1e-8

        elif args["end_chr"] == 1 or args["end_chr"] == 2:
            sequence_length = 50_000_000
            rate = 1e-8

        elif args["end_chr"] == 22:
            rate = msprime.RateMap.read_hapmap(args["hapmap_chr1"].replace("chr1", f"chr{chrom}"))
            sequence_length = int(rate.position[-1])

        return sequence_length, rate



class DemographicSetup:
    @staticmethod
    def create(args):

        print("Getting custom demographic model")
        spec = importlib.util.spec_from_file_location("custom_model", args["custom_demo"]["path"])
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        return getattr(module, args["custom_demo"]["object"])


class Simulation:
    @staticmethod
    def create(args, chrom, path):

        seed = args["seed"]
        iter_n = args["iter_n"]

        if args["custom_sim"]["path"]:
            print("Custom simulation provided.")
            spec = importlib.util.spec_from_file_location("custom_model", args["custom_sim"]["path"])
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)

            simulate = getattr(module, args["custom_sim"]["object"])

            ts, rate = simulate(chrom, args)
            demography = []

        else:
            sequence_length, rate = GenomeSetup.create(args, chrom)
            demography = DemographicSetup.create(args)

            print(demography.debug())

            if args["pedigree"]["pedigree_mode"]:
                pedigree = msprime.parse_pedigree(open(f"{path}_WF.pedigree"), sequence_length=sequence_length)

                init_state = msprime.sim_ancestry(initial_state=pedigree,
                                model="fixed_pedigree",
                                recombination_rate=rate,
                                end_time=args["pedigree"]["gen_end"],
                                discrete_genome=True,
                                random_seed=seed
                                )
                seed += 10
                init_samples = None

            else:
                init_state = None
                init_samples = args["samples"]

            ts = msprime.sim_ancestry(samples={"pop_0": init_samples} if init_samples else None,
                              initial_state=init_state,
                              sequence_length=sequence_length,
                              discrete_genome=True,
                              recombination_rate=rate,
                              demography=demography,
                              random_seed=seed
                              )
            seed += 20

        return msprime.sim_mutations(ts, rate=1e-8, random_seed=seed - 10), rate, [demography], seed

def run_hapibd(prefix, gb, hapibd_jar='/users/cwilli50/hap-ibd.jar'):
    command = [
        'java',
        '-jar',
        f'-Xmx{int(gb*0.8)}g',
        hapibd_jar,
        f'gt={prefix}.vcf.gz',
        f'map={prefix}.map',
        f'out={prefix}'
    ]
    
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print(f"Command completed successfully: {result.stdout}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e.stderr}")
        return False

def get_tree_index(ts, interval_size):
    breakpoints = ts.breakpoints(as_array=True)
    tree_index = []
    index = None
    cur_interval = interval_size
    while True:
        index = np.argmin(breakpoints < cur_interval) - 1
        cur_interval += interval_size
        if index > -1:
            tree_index.append(index)
        else:
            break
    return tree_index

def get_node(vcf_id, haplotype):
    return int(vcf_id.split("_")[1])*2 + (haplotype==2)

def assign_rca(pair_df, small_segments, large_segments):
    def overlap(x1, x2, y1, y2):
        o = (min(x2, y2) - max(x1, y1))
        return o / (x2 - x1) if o > 0 else 0

    segments = [[seg.left, seg.right, seg.node] for seg in large_segments]
    segments += [[seg.left, seg.right, seg.node] for seg in small_segments if (seg.right-seg.left) < 500_000]
    segments.sort()

    out_mrcas = dict()

    for index, row in pair_df.iterrows():
        nodes = dict()
        for (y1, y2, node) in segments:
            p = overlap(row[5], row[6], y1, y2)
            if p > 1 or p < 0:
                print(row[5], row[6], y1, y2, node)
            nodes[node] = nodes.get(node, 0) + p

        out_mrcas[index] = nodes

    return out_mrcas

class TMRCA:
    def __init__(self):
        self.mrca = dict()
        self.tmrca = dict()
        self.proportion = dict()
        self.node_to_indv = dict()
        self.indv_to_name = dict()

    def add_segment(self, index, mrcas, tmrcas, proportions):
        self.mrca[index] = mrcas
        self.tmrca[index] = tmrcas
        self.proportion[index] = proportions

    def get_info(self, index, attr):
        return getattr(self, attr)[index]

    def add_individuals(self, ts):
        self.node_to_indv = {i: ts.node(i).individual for i in set(it.chain(*self.mrca.values()))}
        try:
            self.indv_to_name = {i: ts.individual(i).metadata["file_id"] if i != -1 else -1 for i in self.node_to_indv.values()}
        except:
            try:
                self.indv_to_name = {i: ts.individual(i).metadata["individual_name"] if i != -1 else -1 for i in self.node_to_indv.values()}
            except:
                self.indv_to_name = dict()

    def save_it(self, prefix):
        with open(f"{prefix}.tmrca.pkl", "wb") as pklf:
            pkl.dump(self, pklf)
        
def add_tmrca(prefix, ts, dump, interval_size=1000000):

    if dump:
        ts.dump(prefix + ".trees")

    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=FutureWarning)
        df = pd.read_csv(prefix + ".ibd.gz", sep="\\s+", header=None)

    df["id1"] = df[[0, 1]].apply(lambda x: get_node(*x), axis=1)
    df["id2"] = df[[2, 3]].apply(lambda x: get_node(*x), axis=1)
    df["phy_len"] = df[6] - df[5]

    small_segments = df[df.phy_len < 500_000]
    large_segments = df[df.phy_len >= 500_000]

    small_ibd_obj = ts.ibd_segments(within=[int(i) for i in set(it.chain(*small_segments[["id1","id2"]].values))],
                                    min_span=small_segments.phy_len.min()-1,
                                    store_segments=True)

    large_ibd_obj = ts.ibd_segments(within=[int(i) for i in set(it.chain(*large_segments[["id1","id2"]].values))],
                                    min_span=500_000-1,
                                    store_segments=True)

    tmrca = TMRCA()

    for pair, pair_df in df.groupby(["id1", "id2"]):

        mrcas = assign_rca(pair_df, small_ibd_obj.get(pair, []), large_ibd_obj.get(pair, []))

        for index, mrca in mrcas.items():
            rcas = list(mrca.keys())

            tmrca.add_segment(index,
                              mrcas=rcas,
                              tmrcas=list(ts.nodes_time[rcas]),
                              proportions=[mrca[i] for i in rcas])

    tmrca.add_individuals(ts)

    tmrca.save_it(prefix)


def base_seed(path, iter_n):

    hash_object = hashlib.md5(path.encode())

    full_hex = hash_object.hexdigest()

    return (int(full_hex[-16:], 16) + iter_n)  % (2**32)

def sim(path, iter_n, chrom):

    iteration_seed = base_seed(path, int(iter_n))
    print(f"Iteration seed: {iteration_seed}")

    try:
        seed = int(sys.argv[4])
    except:
        seed = np.random.randint(1, 1000000)

    yargs = yaml.safe_load(open(f"{path}/args.yaml"))

    config = load_config()

    yargs["seed"] = seed
    yargs["iteration_seed"] = iteration_seed
    yargs["iter_n"] = iter_n
    yargs["chrom"] = chrom
 
    print(f"Random seed: {seed}")


    # Dead code for direct CLI use only
    if int(chrom) == 0 and yargs["pedigree"]["pedigree_mode"]:
        log = open(f"{path}/simulation.log", "w")

        demography = DemographicSetup.create(yargs)

        gen_arr, Ne_arr = WFSetup.create(yargs, path, demography)

        print("Done creating WF pedigree file.")

        for g, Ne in zip(gen_arr, Ne_arr):
            log.write(f"time {g}: {Ne} individuals\n")

        log.close()

        sys.exit()

    # if random.random() < 0.25: raise RuntimeError("Random crash!")

    ts, rate, demography, seed = Simulation.create(yargs, chrom, f"{path}/iter{iter_n}")


    # log.write(f"\n\n{demography[0].debug() if len(demography)>0 else 'ts supplied!'}\n")


    prefix = f"{path}/iter{iter_n}_chr{chrom}"

    write_vcf(ts, prefix, chrom, rate, seed, snps_pkl=config["maf_pickle"])

    run_hapibd(prefix, yargs["gb"], hapibd_jar=config["hap_ibd_jar"])

    add_tmrca(prefix, ts, False)

    if os.path.exists(f"{prefix}.ibd.gz"):

        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=FutureWarning)

            df = pd.read_csv(f"{prefix}.ibd.gz", nrows=10, sep="\\s+", header=None)

        if df.shape[0] == 10:

            os.remove(f"{prefix}.vcf.gz")
            os.remove(f"{prefix}.hbd.gz")

            print("Success!")

            sys.exit()

    else:
        err = open(f"{path}/errors/iter{iter_n}_chr{chrom}.err", "w")
        err.write(f"No .ibd.gz file found for iter {iter_n}, chromosome {chrom}")
        err.close()

    return True

if __name__ == "__main__":

    path = sys.argv[1]
    iter_n = sys.argv[2]
    chrom = sys.argv[3]

    sim(path, iter_n, chrom)

    

