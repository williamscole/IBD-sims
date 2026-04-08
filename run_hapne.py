import os
import shutil
import tempfile
import numpy as np
from configparser import ConfigParser
import sys
import pandas_plink._read as _plink_read
import pandas as pd

_original_read_csv = pd.read_csv

# Monkey patch: delim_whitespace is deprecated in newer versions
def _patched_read_csv(*args, **kwargs):
    if kwargs.pop("delim_whitespace", False):
        kwargs["sep"] = "\\s+"
    return _original_read_csv(*args, **kwargs)

pd.read_csv = _patched_read_csv
_plink_read.read_csv = _patched_read_csv

import hapne.convert.tools
import hapne.ld as hapne_ld_module
from hapne import hapne_ld as _hapne_ld


def hapne_tmp_map(input_map: str) -> tuple[str, str]:
    """Convert a plink .map file to per-chromosome SHAPEIT-format map files.

    Plink .map columns: chrom, snp_id, genetic_pos (cM), physical_pos (bp)
    SHAPEIT columns:    physical_pos (bp), rate (cM/Mb), genetic_pos (cM)

    Returns (tmp_dir, pattern) where pattern is like '/tmp/xyz/chr@.txt'.
    The caller is responsible for cleaning up tmp_dir via shutil.rmtree.
    """
    df = pd.read_csv(
        input_map,
        sep="\\s+",
        header=None,
        names=["chrom", "snp_id", "cm", "bp"],
    )

    tmp_dir = tempfile.mkdtemp()

    for chrom, chrom_df in df.groupby("chrom"):
        chrom_df = chrom_df.sort_values("bp").reset_index(drop=True)

        bp = chrom_df["bp"].to_numpy(dtype=float)
        cm = chrom_df["cm"].to_numpy(dtype=float)

        rate = np.zeros(len(bp))
        delta_bp = np.diff(bp)
        delta_cm = np.diff(cm)

        with np.errstate(invalid="ignore", divide="ignore"):
            rate[1:] = np.where(delta_bp > 0, delta_cm / (delta_bp / 1e6), 0.0)

        out_df = pd.DataFrame({
            "bp": bp.astype(int),
            "rate": rate,
            "cm": cm,
        })

        out_path = os.path.join(tmp_dir, f"chr{chrom}.txt")
        out_df.to_csv(out_path, sep=" ", index=False, header=False)

    pattern = os.path.join(tmp_dir, "chr@.txt")
    return tmp_dir, pattern


def _make_sim_regions(end_chr: int, chr_len_bp: int = 100_000_000) -> pd.DataFrame:
    """Build a HapNe-compatible regions DataFrame for simulated chromosomes.

    Each chromosome is split into two equal arms (p and q).
    Returns a DataFrame with columns: CHR, FROM_BP, TO_BP, NAME.
    """
    half = chr_len_bp // 2
    rows = []
    for chrom in range(1, end_chr + 1):
        rows.append({
            "CHR": chrom,
            "FROM_BP": 1,
            "TO_BP": half,
            "NAME": f"chr{chrom}p",
        })
        rows.append({
            "CHR": chrom,
            "FROM_BP": half + 1,
            "TO_BP": chr_len_bp,
            "NAME": f"chr{chrom}q",
        })
    return pd.DataFrame(rows)


def _make_sim_map_files(
    input_map: str,
    regions: pd.DataFrame,
    output_folder: str,
):
    """Write per-arm SHAPEIT map files into output_folder/DATA/.

    Files are named to match HapNe's fallback convention:
    chr{N}.from{from_bp}.to{to_bp}.map
    """
    data_dir = os.path.join(output_folder, "DATA")
    os.makedirs(data_dir, exist_ok=True)

    df = pd.read_csv(
        input_map,
        sep="\\s+",
        header=None,
        names=["chrom", "snp_id", "cm", "bp"],
    )

    for _, region in regions.iterrows():
        chrom = region["CHR"]
        from_bp = region["FROM_BP"]
        to_bp = region["TO_BP"]

        arm_df = df[
            (df["chrom"] == chrom) &
            (df["bp"] >= from_bp) &
            (df["bp"] <= to_bp)
        ].sort_values("bp").reset_index(drop=True)

        bp = arm_df["bp"].to_numpy(dtype=float)
        cm = arm_df["cm"].to_numpy(dtype=float)

        rate = np.zeros(len(bp))
        if len(bp) > 1:
            delta_bp = np.diff(bp)
            delta_cm = np.diff(cm)
            with np.errstate(invalid="ignore", divide="ignore"):
                rate[1:] = np.where(delta_bp > 0, delta_cm / (delta_bp / 1e6), 0.0)

        out_df = pd.DataFrame({"bp": bp.astype(int), "rate": rate, "cm": cm})
        out_path = os.path.join(data_dir, f"chr{chrom}.from{from_bp}.to{to_bp}.map")
        out_df.to_csv(out_path, sep=" ", index=False, header=False)


def run_hapne_ld(
    vcf_file: str,
    input_map: str,
    output_folder: str,
    population_name: str,
    end_chr: int,
    chr_len_bp: int = 100_000_000,
    genome_build: str = "grch37",
):
    """Run the full HapNe-LD pipeline on simulated data."""

    if end_chr != 22:
        regions = _make_sim_regions(end_chr, chr_len_bp)
        _make_sim_map_files(input_map, regions, output_folder)

        def _get_regions(build):
            return regions

        def _get_region(index, build):
            return regions.iloc[index]

        hapne.convert.tools.get_regions = _get_regions
        hapne.convert.tools.get_region = _get_region
        hapne_ld_module.get_regions = _get_regions
        hapne_ld_module.get_region = _get_region

        # map config value is ignored by HapNe for the simulated case —
        # it builds the path itself from output_folder/DATA/chr@.from.to.map
        map_pattern = os.path.join(output_folder, "DATA", "chr@.map")
        tmp_dir = None

    else:
        tmp_dir, map_pattern = hapne_tmp_map(input_map)

    try:
        hapne_config = ConfigParser()
        hapne_config["CONFIG"] = {
            "vcf_file": vcf_file,
            "map": map_pattern,
            "output_folder": output_folder,
            "population_name": population_name,
            "genome_build": genome_build,
            "pseudo_diploid": "False",
        }

        print("[HapNe-LD] stage 1: split/convert VCF")
        hapne.convert.tools.split_convert_vcf(hapne_config)

        print("[HapNe-LD] stage 2: compute LD")
        hapne_ld_module.compute_ld(hapne_config)
        hapne_ld_module.compute_ccld(hapne_config)

        print("[HapNe-LD] stage 3: HapNe-LD")
        _hapne_ld(hapne_config)

    finally:
        if tmp_dir is not None:
            shutil.rmtree(tmp_dir)

if __name__ == "__main__":

    run_hapne_ld(vcf_file=sys.argv[1],
                 input_map=sys.argv[2],
                 output_folder=sys.argv[3],
                 population_name=sys.argv[4],
                 end_chr=int(sys.argv[5]))