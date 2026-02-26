import pandas as pd
import numpy as np
import pickle as pkl
import argparse

BUCKETS = [(-1, 0), (0, 0.01), (0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5)]


def maf_bucket(chrom, afreq_chr1):
    afreq_path = afreq_chr1.replace("chr1", f"chr{chrom}")
    counts = np.zeros(len(BUCKETS))
    df = pd.read_csv(afreq_path, delim_whitespace=True)
    for index, (i, j) in enumerate(BUCKETS):
        counts[index] = df.MAF.apply(lambda x: i < x <= j).sum()
    return counts


def get_snp_density(chrom, bim_chr1, w=1_000_000, n=100):
    bim_path = bim_chr1.replace("chr1", f"chr{chrom}")
    bim = pd.read_csv(bim_path, delim_whitespace=True, header=None)
    max_end = bim.iloc[-1][3] - w
    min_end = bim.iloc[0][3] + 1
    random_windows = np.random.randint(min_end, max_end, n)
    densities = []
    for w1 in random_windows:
        n_snps = bim[3].apply(lambda x: w1 <= x < (w1 + w)).sum()
        densities.append(n_snps / w)
    return densities


class SNPs:
    def __init__(self, bucket_counts, densities):
        self.freqs = np.array(bucket_counts) / sum(bucket_counts)
        self.buckets = BUCKETS
        self.n_buckets = len(BUCKETS)
        self.density = np.mean(densities)

    def get_counts(self, chrom_len):
        n_snps = int(chrom_len * self.density)
        buckets = []
        for index in range(self.n_buckets):
            freq = self.freqs[index]
            bucket_count = int(freq * n_snps)
            buckets.append([*self.buckets[index], bucket_count])
        return buckets

    def save_me(self, output):
        with open(output, "wb") as pklf:
            pkl.dump(self, pklf)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Build a SNP density/MAF distribution pickle from plink allele frequency "
                    "and bim files. Expects one file per chromosome; provide the chr1 path and "
                    "all others are inferred by replacing 'chr1' with 'chr{n}'. "
                    "afreq files should be plink2 --freq output; "
                    "bim files should be standard plink .bim format."
    )
    parser.add_argument("--afreq-chr1", type=str, required=True,
                        help="Path to chr1 allele frequency file (plink2 --freq output).")
    parser.add_argument("--bim-chr1", type=str, required=True,
                        help="Path to chr1 bim file (plink .bim format).")
    parser.add_argument("--output", type=str, required=True,
                        help="Path to write the output pickle.")
    args = parser.parse_args()

    counts = np.zeros(len(BUCKETS))
    snp_densities = []

    for chrom in range(1, 23):
        new_counts = maf_bucket(chrom, args.afreq_chr1)
        counts += new_counts
        snp_densities.extend(get_snp_density(chrom, args.bim_chr1))

    snp_counts = SNPs(counts, snp_densities)
    snp_counts.save_me(args.output)
    print(f"Wrote: {args.output}")