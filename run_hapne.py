import os
import shutil
import tempfile
import numpy as np
from configparser import ConfigParser
from concurrent.futures import ProcessPoolExecutor
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


# ── Region helpers ────────────────────────────────────────────────────────────

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
    """Write per-arm SHAPEIT map files into output_folder/DATA/."""
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


# ── Map conversion ────────────────────────────────────────────────────────────

def hapne_tmp_map(input_map: str) -> tuple[str, str]:
    """Convert a plink .map file to per-chromosome SHAPEIT-format map files."""
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


# ── Parallel worker functions ─────────────────────────────────────────────────
# Each worker runs in a separate process, so monkey-patches and config
# must be re-applied. We pass serialisable arguments (config dict, regions
# as a list of dicts) rather than ConfigParser objects.

def _apply_patches(regions_dicts=None):
    """Re-apply monkey-patches in a worker process."""
    import pandas as pd
    import pandas_plink._read as _plink_read

    _orig = pd.read_csv.__wrapped__ if hasattr(pd.read_csv, "__wrapped__") else None

    # Only patch if not already patched
    if _orig is None:
        _orig = pd.read_csv

        def _patched(*args, **kwargs):
            if kwargs.pop("delim_whitespace", False):
                kwargs["sep"] = "\\s+"
            return _orig(*args, **kwargs)

        _patched.__wrapped__ = _orig  # mark as patched
        pd.read_csv = _patched
        _plink_read.read_csv = _patched

    if regions_dicts is not None:
        regions_df = pd.DataFrame(regions_dicts)

        def _get_regions(build="grch37"):
            return regions_df

        def _get_region(index, build="grch37"):
            return regions_df.iloc[index]

        import hapne.convert.tools
        import hapne.ld as _ld

        hapne.convert.tools.get_regions = _get_regions
        hapne.convert.tools.get_region = _get_region
        _ld.get_regions = _get_regions
        _ld.get_region = _get_region


def _rebuild_config(config_dict: dict) -> ConfigParser:
    """Rebuild a ConfigParser from a plain dict."""
    config = ConfigParser()
    config["CONFIG"] = config_dict
    return config


def _worker_compute_ld(args):
    """Worker: compute LD for one region."""
    region_index, config_dict, regions_dicts = args
    _apply_patches(regions_dicts)
    config = _rebuild_config(config_dict)
    hapne_ld_module.compute_ld_in_parallel(region_index, config)
    return region_index


def _worker_compute_ccld(args):
    """Worker: compute CCLD for one pair, write result to a temp file."""
    job_index, config_dict, regions_dicts, tmp_dir = args
    _apply_patches(regions_dicts)
    config = _rebuild_config(config_dict)

    # Compute the CCLD quantities (reimplemented to avoid file-append race)
    build = config.get('CONFIG', 'genome_build', fallback='grch37')
    reg1, reg2 = hapne_ld_module.get_regions_from_index(job_index, build)
    maf = config.getfloat("CONFIG", "maf", fallback=0.25)
    pseudo_diploid = config.getboolean("CONFIG", "pseudo_diploid", fallback=False)

    folder = hapne_ld_module.get_genotypes_location(config)
    genotype1, _ = hapne_ld_module.load_and_preprocess_file(folder + "/" + reg1 + ".bed", maf)
    genotype2, _ = hapne_ld_module.load_and_preprocess_file(folder + "/" + reg2 + ".bed", maf)

    ccld, expected_ccld, bessel_factor, s_corr = hapne_ld_module._compute_cc_quantities(
        genotype1, genotype2, int(1e6), pseudo_diploid
    )

    # Write to individual temp file instead of appending to shared file
    out_path = os.path.join(tmp_dir, f"ccld_{job_index}.csv")
    with open(out_path, "w") as f:
        f.write(f"{reg1},{reg2},{ccld},{expected_ccld},{bessel_factor},{s_corr}\n")

    return job_index


# ── Main pipeline ─────────────────────────────────────────────────────────────

def run_hapne_ld(
    vcf_file: str,
    input_map: str,
    output_folder: str,
    population_name: str,
    end_chr: int,
    chr_len_bp: int = 100_000_000,
    genome_build: str = "grch37",
    workers: int = 1,
):
    """Run the full HapNe-LD pipeline on simulated data."""

    regions_dicts = None  # only set for non-standard (simulated) chromosomes

    if end_chr != 22:
        regions = _make_sim_regions(end_chr, chr_len_bp)
        _make_sim_map_files(input_map, regions, output_folder)
        regions_dicts = regions.to_dict("records")

        def _get_regions(build="grch37"):
            return regions

        def _get_region(index, build="grch37"):
            return regions.iloc[index]

        hapne.convert.tools.get_regions = _get_regions
        hapne.convert.tools.get_region = _get_region
        hapne_ld_module.get_regions = _get_regions
        hapne_ld_module.get_region = _get_region

        map_pattern = os.path.join(output_folder, "DATA", "chr@.map")
        tmp_dir_map = None

    else:
        tmp_dir_map, map_pattern = hapne_tmp_map(input_map)

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

        # Serialisable copy for worker processes
        config_dict = dict(hapne_config["CONFIG"])

        if regions_dicts is not None:
            nb_regions = len(regions_dicts)
        else:
            from hapne.utils import get_regions
            nb_regions = get_regions(genome_build).shape[0]

        nb_pairs = nb_regions * (nb_regions - 1) // 2

        # ── Stage 1: split/convert VCF ────────────────────────────────────
        print("[HapNe-LD] stage 1: split/convert VCF")
        hapne.convert.tools.split_convert_vcf(hapne_config)

        if workers <= 1:
            # Sequential (original behaviour)
            print(f"[HapNe-LD] stage 2a: compute LD ({nb_regions} regions, sequential)")
            hapne_ld_module.compute_ld(hapne_config)
            print(f"[HapNe-LD] stage 2b: compute CCLD ({nb_pairs} pairs, sequential)")
            hapne_ld_module.compute_ccld(hapne_config)

        else:
            # ── Stage 2a: LD in parallel ──────────────────────────────────
            print(f"[HapNe-LD] stage 2a: compute LD ({nb_regions} regions, {workers} workers)")
            ld_args = [
                (ii, config_dict, regions_dicts)
                for ii in range(nb_regions)
            ]
            with ProcessPoolExecutor(max_workers=workers) as pool:
                for idx in pool.map(_worker_compute_ld, ld_args):
                    print(f"  LD region {idx + 1}/{nb_regions} done")

            # ── Stage 2b: CCLD in parallel ────────────────────────────────
            # Each worker writes to its own temp file; we concatenate after.
            print(f"[HapNe-LD] stage 2b: compute CCLD ({nb_pairs} pairs, {workers} workers)")
            ccld_tmp_dir = tempfile.mkdtemp()
            try:
                ccld_args = [
                    (ii, config_dict, regions_dicts, ccld_tmp_dir)
                    for ii in range(nb_pairs)
                ]
                with ProcessPoolExecutor(max_workers=workers) as pool:
                    for idx in pool.map(_worker_compute_ccld, ccld_args):
                        if (idx + 1) % 50 == 0 or idx + 1 == nb_pairs:
                            print(f"  CCLD pair {idx + 1}/{nb_pairs} done")

                # Concatenate CCLD temp files into the final .ccld file
                output_ld_folder = hapne_ld_module.get_ld_output_folder(hapne_config)
                analysis_name = hapne_ld_module.get_analysis_name(hapne_config)
                ccld_file = output_ld_folder / f"{analysis_name}.ccld"
                with open(ccld_file, "w") as out:
                    out.write("REGION1, REGION2, CCLD, CCLD_H0, BESSEL_FACTOR, S_CORR\n")
                    for ii in range(nb_pairs):
                        tmp_path = os.path.join(ccld_tmp_dir, f"ccld_{ii}.csv")
                        with open(tmp_path) as f:
                            out.write(f.read())
            finally:
                shutil.rmtree(ccld_tmp_dir)

        # ── Stage 3: fit ──────────────────────────────────────────────────
        print("[HapNe-LD] stage 3: HapNe-LD")
        _hapne_ld(hapne_config)

    finally:
        if tmp_dir_map is not None:
            shutil.rmtree(tmp_dir_map)


if __name__ == "__main__":
    workers = int(sys.argv[6]) if len(sys.argv) > 6 else 1

    run_hapne_ld(
        vcf_file=sys.argv[1],
        input_map=sys.argv[2],
        output_folder=sys.argv[3],
        population_name=sys.argv[4],
        end_chr=int(sys.argv[5]),
        workers=workers,
    )