import pandas as pd
import numpy as np
import pickle as pkl
import sys
import msprime
import itertools as it
from simulations import GenomeSetup

CHROM22 = "/gpfs/data/sramacha/datasets/quebec_sims/simulated_chrom_22.ts.tsz"

PATH = "quebec_trees"

def get_node_table(ts):

    # get the nodes corresponding to the 0th generation and create a pandas dataframe
    node_table = ts.dump_tables().nodes
    nodes = pd.DataFrame({"time": node_table.time, "individual": node_table.individual})

    # get nodes from cur gen
    cur_gen_nodes = nodes[nodes.time==0]

    # get the individual name
    cur_gen_nodes["id1"] = cur_gen_nodes["individual"].apply(lambda x: ts.individual(x).metadata["individual_name"])

    return cur_gen_nodes

def load_chr21():

    import tszip

    print("Loading chromosome 21")
    ts = tszip.decompress(CHROM22.replace("_22", "_21"))
    print("\tLoaded.")

    cur_gen_nodes = get_node_table(ts)

    # get the individuals
    keep_nodes = cur_gen_nodes.drop_duplicates("id1").sample(50000)["id1"]

    # write out
    keep_nodes.to_csv(f"{PATH}/random_ids.txt", header=False, index=False)

def simplify_quebec(chrom, seed=42):
    import tszip

    keep_nodes = pd.read_csv(f"{PATH}/random_ids.txt", header=None)
    keep_ids = keep_nodes[0].apply(str).values 

    print(f"Loading chromosome {chrom}")
    ts = tszip.decompress(CHROM22.replace("_22.", f"_{chrom}."))
    print("\tLoaded.")

    node_df = get_node_table(ts)

    id_to_node = node_df.reset_index().groupby('id1')['index'].min().to_dict()
    chr_nodes_to_keep = list(it.chain(*[[id_to_node[pid], id_to_node[pid]+1] for pid in keep_ids]))

    ts = ts.simplify(samples=chr_nodes_to_keep)

    # Reorder individual table so current-gen individuals are in keep_ids order.
    # ts.simplify preserves original individual table order, which differs across
    # chromosomes — this makes tsk_N VCF headers consistent for bcftools concat.
    tables = ts.dump_tables()
    old_individuals = tables.individuals.copy()

    name_to_idx = {
        str(ts.individual(i).metadata["individual_name"]): i
        for i in range(ts.num_individuals)
        if "individual_name" in ts.individual(i).metadata
    }

    keep_set = set(keep_ids)
    new_order = (
        [name_to_idx[pid] for pid in keep_ids] +
        [i for i in range(ts.num_individuals)
         if str(ts.individual(i).metadata.get("individual_name", "")) not in keep_set]
    )

    old_to_new = np.empty(ts.num_individuals, dtype=np.int32)
    for new_idx, old_idx in enumerate(new_order):
        old_to_new[old_idx] = new_idx

    tables.individuals.clear()
    for old_idx in new_order:
        row = old_individuals[old_idx]
        new_parents = [old_to_new[p] if p >= 0 else -1 for p in row.parents]
        tables.individuals.append(
            flags=row.flags,
            location=row.location,
            parents=new_parents,
            metadata=row.metadata,
        )

    indiv_col = tables.nodes.individual.copy()
    mask = indiv_col >= 0
    indiv_col[mask] = old_to_new[indiv_col[mask]]
    tables.nodes.individual = indiv_col

    ts = tables.tree_sequence()

    tszip.compress(ts, f"{PATH}/chr{chrom}_random_50000.trees.tsz")

def load_random_50000(chrom, args):

    import tszip

    iteration_seed = args.get("iteration_seed", 1000)

    np.random.seed(iteration_seed)

    ts = tszip.decompress(f"{PATH}/chr{chrom}_random_50000.trees.tsz")

    n_samples = ts.num_samples // 2

    expected = 50000
    if n_samples != expected:
        raise ValueError(
            f"chr{chrom} simplified tree has {n_samples} individuals, expected {expected}. "
            f"Re-run simplify_quebec to fix cross-chromosome inconsistency."
        )

    random_nodes = np.random.choice(np.arange(n_samples), args.get("samples", 1000), replace=False) * 2
    random_nodes = np.sort(np.concatenate((random_nodes, random_nodes+1)))

    if len(random_nodes) < (n_samples * 2):
        ts = ts.simplify(samples=random_nodes)

    _, rate = GenomeSetup.create(args, chrom)

    return ts, rate

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

    if sys.argv[1] == "init":
        load_chr21()

    else:
        simplify_quebec(int(sys.argv[1]))


