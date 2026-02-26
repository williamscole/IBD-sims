import numpy as np

from purple import *

if __name__ == "__main__":

    filen="/gpfs/data/sramacha/cwilli50/chapter3/uk_biobank/chr1.ibd.gz"
    # filen = "test/chr1.ibd.gz"

    n_chrom = 22

    file_to_write="uk_biobank/iter1.npy"
    # file_to_write="test/iter1.npy"

    ids_to_keep = np.loadtxt("uk_biobank/trio_ids.txt", dtype=str)
    # ids_to_keep = [f"tsk_{i}" for i in range(0, 1000)]

    nthreads=22
    # nthreads=1

    mat = readin_ibd(filen, n_chrom, ids_to_keep, file_to_write, nthreads=nthreads)