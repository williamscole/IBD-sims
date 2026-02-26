import pandas as pd
import numpy as np
import itertools as it
import yaml
import os
import sys
import importlib


from concat_tmrca import load_tmrca
from purple import readin_ibd




def load_ibd(path, iter_n):

    ibd_df = pd.read_csv(f"{path}/iter{iter_n}.ibd.gz", delim_whitespace=True, header=None).drop_duplicates()

    tmrca_df = load_tmrca(f"{path}/iter{iter_n}.tmrca.gz")

    assert ibd_df.shape[0] == tmrca_df.shape[0]

    assert (tmrca_df["chromosome"].values == ibd_df[4].values).mean() == 1

    return ibd_df.merge(tmrca_df[["proportion","tmrca","mrca","chrom_index"]], left_index=True, right_index=True)


##### For purple nodes

def bin_purple(purple, path=None, iter_n=None):

    if type(path) == str:
        purple = np.load(f"{path}/iter{iter_n}.npy")

    regions = [tuple([i]) for i in range(2, 11)]
    regions += [tuple(range(i, i+5)) for i in range(11, 31, 5)]
    regions += [tuple(range(i, i+10)) for i in range(31, 101, 10)]

    new_purple = np.zeros((len(regions), len(regions), 2))

    for i, j in it.combinations_with_replacement(regions, r=2):
        tot, purp = 0, 0
        for idx1, idx2 in it.product(i, j):
            purp += purple[idx1, idx2, 0]
            tot += purple[idx1, idx2, 1]
        idx = regions.index(i)
        jdx = regions.index(j)
        new_purple[idx, jdx, 0] = purp
        new_purple[idx, jdx, 1] = tot

    return new_purple, [f"({min(i)}, {max(i)})" if len(i) > 1 else f"({i[0]},)" for i in regions]

def stack_purple(path, bin_l=True):

    n_iter = yaml_arg(path, "iter")

    mats = []
    if bin_l:
        for i in range(1, n_iter+1):
            try:
                binned = bin_purple(None, path=f"{path}", iter_n=i)
            except:
                print("No .npy file for iteration " + str(i))
                continue
            mats.append(binned[0])
            labels = binned[1]
    else:
        labels = None
        for i in range(1, n_iter+1):
            try:
                mat = np.load(f"{path}/iter{i}.npy")
            except:
                print("No .npy file for iteration " + str(i))
                continue
            mats.append(mat)

    return np.stack(mats).sum(axis=0), labels

def run_purple(path):

    n_iter = yaml_arg(path, "iter")

    for i in range(1, n_iter+1):

        ibd_file = f"{path}/iter{i}.ibd.gz"

        readin_ibd(ibd_file, file_to_write=ibd_file.replace(".ibd.gz", ".npy"))

        print(f"Done with {ibd_file}")



def yaml_arg(path, *args):

    yargs = yaml.safe_load(open(f"{path}/args.yaml"))

    out = []
    for i in args:
        if type(i) == str:
            out.append(yargs[i])
        elif len(i) == 2:
            out.append(yargs[i[0]][i[1]])
        elif len(i) == 3:
            out.append(yargs[i[0]][i[1]][i[2]])
    return out if len(out) > 1 else out[0]


###### Analysis functions
def tmrca_analysis(path):
    
    n_iter = yaml_arg(path, "iter")

    segment_counts = np.zeros((3,))

    counts = {i: list() for i in range(2, 290)}

    for i in range(1, n_iter+1):

        ibd_df = load_ibd(path, i)

        segment_counts[0] += ibd_df[ibd_df.tmrca.apply(len)==0].shape[0]
        segment_counts[1] += ibd_df[ibd_df.tmrca.apply(len)>1].shape[0]
        segment_counts[2] += ibd_df.shape[0]

        ibd_df = ibd_df[ibd_df.tmrca.apply(len)==1]
        ibd_df["trmca"] = ibd_df.tmrca.apply(lambda x: int(x[0]))

        ls = ibd_df[7].apply(int).values
        tmrcas = ibd_df["trmca"].values

        for l in np.unique(ls):
            counts[l].extend(list(tmrcas[np.where(ls==l)[0]]))

    qs = [0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05]

    df = pd.DataFrame({"length": [i for i,j in counts.items() if len(j)>0]})

    tmrcas = np.array([np.quantile(counts[i], qs) for i in df.length.values])

    for i in range(len(qs)):
        lab = f"{qs[i]*100}%"
        df[lab] = tmrcas[:,i]

    print(segment_counts)

    df.to_csv(f"{path}/tmrca_by_length.csv", sep=",", index=False)

def tmrca_by_kinship(path):

    from kinship import DEGREES

    def degree_of_r(k):
        for d, (k1, k2) in DEGREES.items():
            if k1 <= k < k2:
                return d
        return "11th"

    n_iter = yaml_arg(path, "iter")

    if yaml_arg(path, "end_chr") == 30:
        cm = 3000

    elif yaml_arg(path, "end_chr") == 22:
        cm = 0
        for chrom in range(1, 23):
            tmp = pd.read_csv(f"/users/cwilli50/reference/geneticMap-GRCh37/genetic_map_GRCh37_chr{chrom}.txt.gz", delim_whitespace=True)
            cm += (tmp.iloc[-1]["Map(cM)"] - tmp.iloc[0]["Map(cM)"])

    degree_d = {d: {l: [] for l in range(2, 300)} for d in DEGREES.keys()}

    for i in range(1, n_iter+1):
        print(f"Iteration: {i}")

        ibd_df = load_ibd(path, i)

        ibd_df = ibd_df[ibd_df.tmrca.apply(len)==1][[0, 2, 7, "tmrca"]]

        for _, pair_df in ibd_df.groupby([0, 2]):
            k = pair_df[7].sum() / cm

            d = degree_of_r(k)

            for row in pair_df.itertuples():
                degree_d[d][int(row._3)].append(row.tmrca[0])
            
    qs = [0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05]

    mat = np.zeros((len(degree_d.keys()), 300, len(qs)+1))

    sorted_degree = sorted(degree_d.keys())

    for index, d in enumerate(sorted_degree):
        for l, ts in degree_d[d].items():
            length = len(ts)
            if length > 0:
                mat[index, l, 1:] = np.quantile(ts, qs)
            mat[index, l, 0] = length

    np.save(f"{path}/tmrca_by_kinship.npy", mat)


def get_path(path):
    if not os.path.exists(path):
        tmp = pd.read_csv("focus_simulations.tsv", delim_whitespace=True)
        arr = int(path)
        path = tmp.iloc[arr]["path"]
        print(f"Array job. Running array {arr} for {path}")
    return path

def import_from_string(module_name, object_name):
    # Import the module
    module = importlib.import_module(module_name)
    
    # Get the object from the module
    obj = getattr(module, object_name)
    
    return obj

def get_ne(demo_obj, max_g):
    gen_arr = np.arange(0, max_g+1)

    ne_arr = demo_obj.debug().population_size_trajectory(gen_arr)

    return pd.DataFrame({"GEN": gen_arr, "NE": ne_arr[:,0]})
    
###### Get the most likely TMRCA of IBD segment lengths as a funciton of Ne
def tmrca_ne(demo_name, demo_path="demography", max_g=300, min_l=2, max_l=100, step=1):
    def coalescent_prob(g, Ne_arr):
        log_prod = 0
        for g_star in range(1, g):
            log_prod += np.log(1 - 1/(2*Ne_arr[g_star]))
        
        # Add the final term separately
        log_final = -np.log(2*Ne_arr[g])
        
        return log_prod + log_final

    def tmrca_g(g, l, Ne_arr):
        term1 = 2*np.log(g / 50)
        term2 = (-l * g / 50)
        term3 = coalescent_prob(g, Ne_arr)

        return term1 + term2 + term3

    demo = import_from_string(demo_path, demo_name)

    ne_df = get_ne(demo, max_g)

    l = min_l
    rows = []
    while l <= max_l:
        likelis = [tmrca_g(g, l, ne_df["NE"].values) for g in range(1, max_g+1)]
        most_likely_g = np.argmax(likelis) + 1
        rows.append([l, most_likely_g])
        l += step

    df = pd.DataFrame(rows, columns=["length", "gen"])

    return df



if __name__ == "__main__":
    if sys.argv[1] == "tmrca_age":

        path = get_path(sys.argv[2])

        tmrca_analysis(path)

    if sys.argv[1] == "purple_matrix":

        path = get_path(sys.argv[2])

        run_purple(path)

    if sys.argv[1] == "tmrca_length":

        tmrca_ne(sys.argv[2])

    if sys.argv[1] == "tmrca_kinship":

        tmrca_by_kinship(sys.argv[2])

