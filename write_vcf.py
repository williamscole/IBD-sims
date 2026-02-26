import msprime
import pandas as pd
import gzip
import numpy as np
import pickle as pkl

from maf_buckets import SNPs, BUCKETS

def allele_frequencies(ts, sample_sets=None):
    if sample_sets is None:
       sample_sets = [ts.samples()] 
    n = np.array([len(x) for x in sample_sets])
    def f(x):
       return x / n
    return ts.sample_count_stat(sample_sets, f, len(sample_sets), windows='sites', polarised=True, mode='site', strict=False, span_normalise=False)

def physical_to_cm(positions, rates):
    # Calculate interval sizes
    interval_sizes = np.diff(positions)  # 165357 elements
    
    cm_pos, mb_pos = np.cumsum((interval_sizes * rates)[1:])*100, positions[2:]

    return np.concatenate((np.array([0]), cm_pos)), np.concatenate((np.array([0]), mb_pos))

def write_vcf(ts, output, chrom, rate, seed):
    with open("ukb_snps.pkl", "rb") as pklf:
        snps = pkl.load(pklf)

    # Get allele freqs
    af = allele_frequencies(ts)
    af = np.where(af < 0.5, af, 1 - af)

    # Get the sites
    positions = ts.tables.sites.position.astype(int)
    sites_df = pd.DataFrame({"position": positions})
    all_sites = set(sites_df.index)
    # drop multi allelic sites
    mut_df = pd.DataFrame({"site": ts.tables.mutations.site.astype(int), "time": ts.tables.mutations.time}).drop_duplicates("site", keep=False)

    keep_sites = set()
    # MAF bucket counts
    snps.density *= 2
    bucket_counts = snps.get_counts(positions[-1])
    # Iterate through each AF
    for af1, af2, count in bucket_counts:
        np.random.seed(seed)
        # Get indices in the index
        pot_index = np.where(np.logical_and((af>af1).flatten(), (af<=af2).flatten(), positions>0))[0]
        # Get the sites
        pot_sites = mut_df[mut_df.site.isin(pot_index)].values[:,0]
        # Random subset of the sites
        keep_sites |= set(np.random.choice(pot_sites, min(pot_sites.shape[0], count), replace=False))
        seed += 1

    # pot_index = np.where(np.logical_and((af>0.05).flatten(), (af<0.95).flatten(), positions>0))[0]
    # pot_sites = mut_df[mut_df.site.isin(pot_index)]["site"].values

    # np.random.seed(seed)

    # keep_index = np.random.choice([True, False], len(pot_sites), p=[0.10, 0.90])
    # keep_sites = set(pot_sites[keep_index])

    # get the sites to delete
    delete_sites = list(set(all_sites-keep_sites))

    out_ts = ts.delete_sites(delete_sites)

    map_df = pd.DataFrame({"bp": out_ts.sites_position.astype(int)})
    map_df["chrom"] = 1
    map_df["rsid"] = map_df.bp.apply(lambda x: f"chr{chrom}_{x}")

    if rate == 1e-8:
        map_df["cm"] = map_df.bp / 1_000_000

    else:
        cm_pos, mb_pos = physical_to_cm(rate.position, rate.rate)
        map_df["cm"] = np.interp(map_df["bp"], mb_pos, cm_pos)
 
    map_df[["chrom", "rsid", "cm", "bp"]].to_csv(f"{output}.map", header=False, index=False, sep=" ")

    with gzip.open(f"{output}.vcf.gz", "wt") as f:
        out_ts.write_vcf(f)

    print(f"Wrote VCF: {output}.vcf.gz")
