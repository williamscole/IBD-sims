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
import hapne.ibd as hapne_ibd_module
import hapne.backend.IO as hapne_io_module
from hapne import hapne_ld as _hapne_ld
from hapne import hapne_ibd as _hapne_ibd


# ── Region helpers ────────────────────────────────────────────────────────────

def _make_sim_regions(end_chr: int, input_map: str) -> pd.DataFrame:
    """Build a HapNe-compatible regions DataFrame for simulated chromosomes.

    Infers each chromosome's length from the map file, then splits into
    two equal arms (p and q).
    Returns a DataFrame with columns: CHR, FROM_BP, TO_BP, NAME, LENGTH.
    LENGTH is the genetic length of the arm in cM (needed by HapNe-IBD).
    """
    map_df = pd.read_csv(
        input_map, sep="\\s+", header=None, names=["chrom", "snp_id", "cm", "bp"],
    )

    rows = []
    for chrom in range(1, end_chr + 1):
        chrom_data = map_df.loc[map_df["chrom"] == chrom]
        if chrom_data.empty:
            raise ValueError(f"Chromosome {chrom} not found in map file: {input_map}")
        max_bp = int(chrom_data["bp"].max())
        half = max_bp // 2

        # Genetic length of each arm in cM
        p_arm = chrom_data[chrom_data["bp"] <= half]
        q_arm = chrom_data[chrom_data["bp"] > half]
        p_len_cm = float(p_arm["cm"].max() - p_arm["cm"].min()) if len(p_arm) > 1 else 0.0
        q_len_cm = float(q_arm["cm"].max() - q_arm["cm"].min()) if len(q_arm) > 1 else 0.0

        rows.append({
            "CHR": chrom,
            "FROM_BP": 1,
            "TO_BP": half,
            "NAME": f"chr{chrom}p",
            "LENGTH": p_len_cm,
        })
        rows.append({
            "CHR": chrom,
            "FROM_BP": half + 1,
            "TO_BP": max_bp,
            "NAME": f"chr{chrom}q",
            "LENGTH": q_len_cm,
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
        import hapne.backend.IO as _io

        hapne.convert.tools.get_regions = _get_regions
        hapne.convert.tools.get_region = _get_region
        _ld.get_regions = _get_regions
        _ld.get_region = _get_region
        _io.get_regions = _get_regions
        _io.get_region = _get_region


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
    genome_build: str = "grch37",
    workers: int = 1,
):
    """Run the full HapNe-LD pipeline on simulated data."""

    regions_dicts = None  # only set for non-standard (simulated) chromosomes

    if end_chr != 22:
        regions = _make_sim_regions(end_chr, input_map)
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
        hapne_io_module.get_regions = _get_regions
        hapne_io_module.get_region = _get_region

        # Override split_convert_vcf: the original hardcodes range(39) and
        # index bounds, and the plink command lacks --chr-set for non-human data.
        def _custom_split_convert_vcf(config: ConfigParser) -> None:
            save_in = config.get("CONFIG", "output_folder") + "/DATA/GENOTYPES"
            os.makedirs(save_in, exist_ok=True)

            vcf_loc = hapne.convert.tools.get_vcf_location(config)
            keep = config.get("CONFIG", "keep", fallback=None)
            keep_command = f"--keep {keep}" if keep else ""

            for ii in range(len(regions)):
                region = regions.iloc[ii]
                chrom = region["CHR"]
                chr_from = region["FROM_BP"]
                chr_to = region["TO_BP"]

                fallback = config.get("CONFIG", "output_folder") + f"/DATA/chr@.from{chr_from}.to{chr_to}.map"
                map_loc = hapne.convert.tools.get_map_location(config, fallback)

                command = (
                    f"plink --vcf {vcf_loc}.vcf.gz "
                    f"--make-bed "
                    f"{keep_command} "
                    f"--out {save_in}/{region['NAME']} "
                    f"--cm-map {map_loc} "
                    f"--chr {chrom} "
                    f"--from-bp {chr_from} "
                    f"--to-bp {chr_to} "
                    f"--const-fid "
                    f"--threads 1 "
                    f"--memory 2048 "
                    f"--maf 0.249 "
                    f"--snps-only "
                    f"--geno 0.8 "
                    f"--chr-set {end_chr} "
                )
                os.system(command)

        hapne.convert.tools.split_convert_vcf = _custom_split_convert_vcf

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
            # Pre-create output directory to avoid race condition between workers
            ld_output = hapne_ld_module.get_ld_output_folder(hapne_config)

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
        print(f"  DEBUG IO.get_region(0) = {hapne_io_module.get_region(0, 'grch37')['NAME']}")  # DEBUG
        print(f"  DEBUG IO.get_region(1) = {hapne_io_module.get_region(1, 'grch37')['NAME']}")  # DEBUG
        _hapne_ld(hapne_config)

    finally:
        if tmp_dir_map is not None:
            shutil.rmtree(tmp_dir_map)


# ── HapNe-IBD pipeline ───────────────────────────────────────────────────────

# Patch:
# sed -i 's/times\[ii + 1\] = t_quantile/times[ii + 1] = np.asarray(t_quantile).item()/' \
#  /users/cwilli50/.conda/envs/ibd-sims/lib/python3.12/site-packages/hapne/backend/DemographicHistory.py

# sed -i 's/n = n.ravel()/n = np.asarray(n).ravel()/' \
#  /users/cwilli50/.conda/envs/ibd-sims/lib/python3.12/site-packages/hapne/utils.py

def _split_ibd_by_arm(ibd_file: str, regions: pd.DataFrame, output_folder: str):
    """Split a whole-genome IBD file into per-arm IBD files.

    Each IBD segment is assigned to the arm containing its midpoint.
    Output files are named {region_name}.ibd.gz in output_folder.

    IBD columns: id1 hap1 id2 hap2 chrom start_bp end_bp length_cM
    """
    import gzip

    ibd_dir = os.path.join(output_folder, "IBD")
    os.makedirs(ibd_dir, exist_ok=True)

    # Build lookup: (chrom) -> list of (from_bp, to_bp, region_name)
    arm_lookup = {}
    for _, region in regions.iterrows():
        chrom = region["CHR"]
        arm_lookup.setdefault(chrom, []).append(
            (region["FROM_BP"], region["TO_BP"], region["NAME"])
        )

    # Open output file handles
    handles = {}
    for _, region in regions.iterrows():
        path = os.path.join(ibd_dir, f"{region['NAME']}.ibd.gz")
        handles[region["NAME"]] = gzip.open(path, "wt")

    try:
        # Read input (may be gzipped)
        if ibd_file.endswith(".gz"):
            f_in = gzip.open(ibd_file, "rt")
        else:
            f_in = open(ibd_file, "r")

        with f_in:
            for line in f_in:
                parts = line.strip().split()
                if len(parts) < 8:
                    continue
                chrom = int(parts[4])
                start_bp = int(parts[5])
                end_bp = int(parts[6])
                midpoint = (start_bp + end_bp) // 2

                arms = arm_lookup.get(chrom, [])
                for from_bp, to_bp, name in arms:
                    if from_bp <= midpoint <= to_bp:
                        # Write as tab-separated (HapNe's awk uses -F"\t")
                        handles[name].write("\t".join(parts) + "\n")
                        break
    finally:
        for h in handles.values():
            h.close()

    return ibd_dir


def run_hapne_ibd(
    ibd_file: str,
    input_map: str,
    output_folder: str,
    population_name: str,
    nb_samples: int,
    end_chr: int,
    genome_build: str = "grch37",
):
    """Run the full HapNe-IBD pipeline on simulated data.

    Parameters
    ----------
    ibd_file : path to the concatenated .ibd.gz file
    input_map : path to the plink .map file
    output_folder : where to write HapNe output
    population_name : label for the analysis
    nb_samples : number of diploid individuals
    end_chr : number of chromosomes (22 = real human, other = simulated)
    genome_build : genome build (only used when end_chr=22)
    """

    if end_chr != 22:
        regions = _make_sim_regions(end_chr, input_map)

        def _get_regions(build="grch37"):
            return regions

        def _get_region(index, build="grch37"):
            return regions.iloc[index]

        hapne_ibd_module.get_regions = _get_regions
        hapne_ibd_module.get_region = _get_region
        hapne_io_module.get_regions = _get_regions
        hapne_io_module.get_region = _get_region

    
    # Get the number of indv in the file
    tmp_ibd_df = pd.read_csv(ibd_file, header=None, sep="\\s+")
    n_samples = len(set(tmp_ibd_df[0].values) | set(tmp_ibd_df[2].values))
    print(f"{n_samples} unique IDs found in IBD file. nb_samples set to {nb_samples}.")

    # Split IBD file by chromosome arm
    print("[HapNe-IBD] stage 1: split IBD by chromosome arm")
    if end_chr != 22:
        nb_regions = len(regions)
    else:
        from hapne.utils import get_regions
        regions = get_regions(genome_build)
        nb_regions = regions.shape[0]

    ibd_dir = _split_ibd_by_arm(ibd_file, regions, output_folder)

    # Build histograms
    print(f"[HapNe-IBD] stage 2: build histograms ({nb_regions} regions)")
    hapne_config = ConfigParser()
    hapne_config["CONFIG"] = {
        "output_folder": output_folder,
        "population_name": population_name,
        "genome_build": genome_build,
        "nb_samples": str(nb_samples),
        "column_cm_length": "8",
        "ibd_files": ibd_dir,
    }

    hapne_ibd_module.build_hist(hapne_config)

    # Fit
    print("[HapNe-IBD] stage 3: HapNe-IBD")
    _hapne_ibd(hapne_config)

if __name__ == "__main__":

    mode = sys.argv[7]

    if mode == "ld":
        print("Running HapNe-LD")

        workers = int(sys.argv[6])

        run_hapne_ld(
            vcf_file=sys.argv[1],
            input_map=sys.argv[2],
            output_folder=sys.argv[3],
            population_name=sys.argv[4],
            end_chr=int(sys.argv[5]),
            workers=workers,
        )
    if mode == "ibd":
        print("Running HapNe-IBD")
        run_hapne_ibd(
            ibd_file=sys.argv[1],
            input_map=sys.argv[2],
            output_folder=sys.argv[3],
            population_name=sys.argv[4],
            nb_samples=int(sys.argv[6]),
            end_chr=int(sys.argv[5])
        )