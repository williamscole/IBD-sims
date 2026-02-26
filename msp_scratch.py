import multiprocessing as mp
import pandas as pd
import numpy as np
import msprime
from functools import partial
import itertools as it

# Simulation parameters
N = 10_000  # Population size
samples = 10000  # Number of samples to draw
recomb_rate = 1e-8  # Human average recombination rate per base pair per generation

chrom_lens = [62_000_000, 46_000_000, 49_000_000]
length = sum(chrom_lens)

position = [0]; rate = []; cur = 0
for i in chrom_lens:
    position.append(cur + i)
    position.append(cur + i + 1)
    rate.append(recomb_rate)
    rate.append(0.5)
    cur += i

rate = rate[:-1]; position = position[:-1]


rate_map = msprime.RateMap(
    position=position,
    rate=rate
)

def write_map(segment_df, prefix):

    tmp = []

    for chrom, chrom_df in segment_df.groupby("chrom"):
        sites = np.array(sorted(list(set(segment_df.start.values) | set(segment_df.end.values)))).astype(int)

        df = pd.DataFrame({"mb": sites, "cm": sites / 1_000_000})

        df["chrom"] = chrom

        df["rsid"] = df.mb.apply(lambda x: f"rs{x}_chr{chrom}")

        tmp.append(df)

    df = pd.concat(tmp)

    df[["chrom", "rsid", "cm", "mb"]].to_csv(f"{prefix}.map", index=False, header=False, sep=" ")

def merge_segments(segments):
    # Sort by start position
    sorted_segs = sorted(segments, key=lambda x: x[0])
    
    if not sorted_segs:
        return []
        
    result = [list(sorted_segs[0])]
    
    for current in sorted_segs[1:]:
        last = result[-1]
        
        # Only check if end of last equals start of current
        if last[1] == current[0]:
            # Keep start from first segment (already there)
            # Update end to second segment's end
            last[1] = current[1]
            # Keep n1, n2 from first segment (already there)
            # Sum the cm values
            last[4] += current[4]
        else:
            # Add as new segment
            result.append(list(current))

    if len(segments) != len(result):
        print("before: ", segments)
        print("after: ", segments)
    
    return result

def get_ibd_segments(ts, samples, chrom):
    ibd_segments = ts.ibd_segments(within=np.arange(samples*2), store_pairs=True, 
                                 store_segments=True, min_span=1_000_000)
    tmp = {}
    for (n1, n2), segs in ibd_segments.items():
        id1 = n1 if n1 % 2 == 0 else (n1 - 1)
        id2 = n2 if n2 % 2 == 0 else (n2 - 1)
        pair = tuple(sorted([id1, id2]))
        if pair not in tmp:
            tmp[pair] = {}
        for seg in segs:
            start, end, n = seg.left, seg.right, seg.node
            cm = (end - start) / 1_000_000
            tmp[pair][n] = tmp[pair].get(n, []) + [[int(start), int(end), n1 % 2, n2 % 2, cm]]

    rows = []
    for (id1, id2), ancestor_segs in tmp.items():
        for _, segs in ancestor_segs.items():
            segs = [[id1, id2, chrom, *i] for i in merge_segments(segs)]
            rows += segs

    df = pd.DataFrame(rows, columns=["id1", "id2", "chrom", "start", "end", "hap1", "hap2", "l"])

    return df[["id1", "hap1", "id2", "hap2", "chrom", "start", "end", "l"]]

def simulate_chromosome(chrom, params):
    """Simulate a single chromosome and get its IBD segments"""
    ts = msprime.sim_ancestry(
        samples=params['samples'],
        sequence_length=params['sequence_length'],
        recombination_rate=params['recomb_rate'],
        population_size=params['Ne'],
        discrete_genome=True,
        end_time=params['end_time']
    )
    return get_ibd_segments(ts, params['samples'], chrom)

def parallel_ajhg_simulation(n_processes=None):
    """Run the AJHG simulation in parallel"""
    if n_processes is None:
        n_processes = mp.cpu_count() - 1  # Leave one core free
    
    # Simulation parameters
    params = {
        'Ne': 10_000,
        'samples': 1_000,
        'recomb_rate': 1e-8,
        'sequence_length': 100_000_000,
        'end_time': 100
    }
    
    # Create a partial function with fixed parameters
    sim_func = partial(simulate_chromosome, params=params)
    
    # Create chromosome numbers (1-30)
    chroms = range(1, 2)
    
    # Run simulations in parallel
    with mp.Pool(processes=n_processes) as pool:
        results = pool.map(sim_func, chroms)
    
    # Combine results
    segments = pd.concat(results)

    write_map(segments, "ajhg_simulation")
    
    # Save results
    segments.to_csv("ajhg_simulation.ibd", header=False, index=False, sep="\t")
    
    return segments




if __name__ == "__main__":
    parallel_ajhg_simulation()
    # # Run the simulation
    # ts = msprime.sim_ancestry(
    #     samples=samples,
    #     sequence_length=length,
    #     recombination_rate=rate_map,
    #     population_size=N,
    #     discrete_genome=True,
    #     end_time=40,
    #     model="dtwf",
    #     random_seed=42  # For reproducibility
    # )

    # # Basic statistics about the simulation
    # print(f"Number of trees: {ts.num_trees}")
    # print(f"Number of samples: {ts.num_samples}")
    # print(f"Sequence length: {ts.sequence_length}")

    # # You can iterate through the trees if needed
    # for tree in ts.trees():
    #     # Just print the first few trees as an example
    #     if tree.index < 3:
    #         print(f"\nTree {tree.index}:")
    #         print(f"Interval: {tree.interval}")
    #         print(f"Total branch length: {tree.total_branch_length}")

    # # Save the tree sequence to a file using command line argument with gzip compression
    # output_filename = sys.argv[-1] + ".trees.gz"
    # ts.dump(output_filename)