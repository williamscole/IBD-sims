import pandas as pd
import os
import sys
import yaml
import numpy as np
import purple
import itertools as it

def main():
    path = sys.argv[-1]

    args = yaml.safe_load(open(f"{path}/args.yaml"))

    n_iter = args["iter"]
    end_chr = args["end_chr"]
    dir_name = args["dir_name"]

    # Get purple node summaries
    mats = []
    for i in range(1, n_iter+1):
        f = f"{path}/{dir_name}/iter{i}.npy"
        if os.path.exists(f):
            mats.append(np.load(f))
        else:
            print(f"File does not exit: {f}")

    comb = np.stack(mats).sum(axis=0)

    purple.plot(comb, label=args["label"], save_to=f"{path}/{dir_name}/purple.png", max_len=99, max_out=True)

    print("Saved purple png to: " + f"{path}/{dir_name}/purple.png")

    lens, tmrcas = [], []
    for i in range(1, n_iter+1):
        f = f"{path}/iter{i}.ibd.gz"
        if os.path.exists(f):
            lens.append(pd.read_csv(f, sep="\\s+", header=None)[7].apply(int).values)
            tmrcas.append(pd.read_csv(f.replace(".ibd.", ".tmrca."), header=None)[0].apply(int).values)
        else:
            print(f"File does not exit: {f}")

    lens = np.array(list(it.chain(*lens)))
    tmrcas = np.array(list(it.chain(*tmrcas)))

    quantile_mat = []
    for i in range(2, np.max(lens)+1):

        tmp = tmrcas[np.where(lens==i)[0]]

        if tmp.shape[0] > 0:
            row = np.quantile(tmp, np.arange(5, 100, 5) / 100)
        else:
            row = 0/np.zeros(19)
        quantile_mat.append(row)

    quantile_mat = np.array(quantile_mat)

    np.save(f"{path}/ibd_segment_ages.npy", quantile_mat)

    print(f"Wrote: {path}/ibd_segment_ages.npy")



if __name__ == "__main__":

    main()