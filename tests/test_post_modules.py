"""
tests/test_post_modules.py

Unit tests for hapne_tmp_map in post_modules.py.
"""

import os
import shutil
import tempfile

import numpy as np
import pandas as pd
import pytest

from post_modules import hapne_tmp_map


# ── Helpers ───────────────────────────────────────────────────────────────────

def write_map(tmp_path, rows):
    """Write a plink .map file from a list of (chrom, snp_id, cm, bp) tuples."""
    path = str(tmp_path / "test.map")
    df = pd.DataFrame(rows, columns=["chrom", "snp_id", "cm", "bp"])
    df.to_csv(path, sep="\t", header=False, index=False)
    return path


def read_shapeit(filepath):
    """Read a SHAPEIT map file back as a DataFrame (bp, rate, cm)."""
    return pd.read_csv(filepath, sep=" ", header=None, names=["bp", "rate", "cm"])


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def single_chr_map(tmp_path):
    """Simple 4-SNP single-chromosome map."""
    rows = [
        (1, "rs1", 0.00, 100_000),
        (1, "rs2", 0.10, 200_000),
        (1, "rs3", 0.25, 350_000),
        (1, "rs4", 0.50, 600_000),
    ]
    return write_map(tmp_path, rows)


@pytest.fixture
def multi_chr_map(tmp_path):
    """3-SNP map across 3 chromosomes."""
    rows = [
        (1, "rs1", 0.00, 100_000),
        (1, "rs2", 0.10, 200_000),
        (1, "rs3", 0.30, 400_000),
        (2, "rs4", 0.00, 150_000),
        (2, "rs5", 0.20, 350_000),
        (3, "rs6", 0.00, 500_000),
        (3, "rs7", 0.05,  600_000),
    ]
    return write_map(tmp_path, rows)


@pytest.fixture
def dup_bp_map(tmp_path):
    """Map with two SNPs sharing the same bp position on one chromosome."""
    rows = [
        (1, "rs1", 0.00, 100_000),
        (1, "rs2", 0.10, 100_000),  # duplicate bp
        (1, "rs3", 0.30, 300_000),
    ]
    return write_map(tmp_path, rows)


@pytest.fixture
def single_snp_map(tmp_path):
    """One SNP per chromosome — minimum possible input."""
    rows = [
        (1, "rs1", 0.00, 100_000),
        (2, "rs2", 0.00, 200_000),
    ]
    return write_map(tmp_path, rows)


# ── Tests ─────────────────────────────────────────────────────────────────────

class TestHapneTmpMapBasic:

    def test_returns_tuple(self, single_chr_map):
        result = hapne_tmp_map(single_chr_map)
        assert isinstance(result, tuple) and len(result) == 2

    def test_pattern_ends_with_wildcard(self, single_chr_map):
        tmp_dir, pattern = hapne_tmp_map(single_chr_map)
        try:
            assert pattern.endswith("chr@.txt")
        finally:
            shutil.rmtree(tmp_dir)

    def test_pattern_is_inside_tmp_dir(self, single_chr_map):
        tmp_dir, pattern = hapne_tmp_map(single_chr_map)
        try:
            assert pattern.startswith(tmp_dir)
        finally:
            shutil.rmtree(tmp_dir)

    def test_output_file_exists(self, single_chr_map):
        tmp_dir, pattern = hapne_tmp_map(single_chr_map)
        try:
            expected = pattern.replace("@", "1")
            assert os.path.exists(expected)
        finally:
            shutil.rmtree(tmp_dir)

    def test_output_has_three_columns(self, single_chr_map):
        tmp_dir, pattern = hapne_tmp_map(single_chr_map)
        try:
            df = read_shapeit(pattern.replace("@", "1"))
            assert df.shape[1] == 3
        finally:
            shutil.rmtree(tmp_dir)

    def test_row_count_matches_input(self, single_chr_map):
        tmp_dir, pattern = hapne_tmp_map(single_chr_map)
        try:
            df = read_shapeit(pattern.replace("@", "1"))
            assert len(df) == 4
        finally:
            shutil.rmtree(tmp_dir)


class TestHapneTmpMapColumns:

    def test_first_snp_rate_is_zero(self, single_chr_map):
        tmp_dir, pattern = hapne_tmp_map(single_chr_map)
        try:
            df = read_shapeit(pattern.replace("@", "1"))
            assert df["rate"].iloc[0] == 0.0
        finally:
            shutil.rmtree(tmp_dir)

    def test_bp_column_matches_input(self, single_chr_map):
        tmp_dir, pattern = hapne_tmp_map(single_chr_map)
        try:
            df = read_shapeit(pattern.replace("@", "1"))
            assert list(df["bp"]) == [100_000, 200_000, 350_000, 600_000]
        finally:
            shutil.rmtree(tmp_dir)

    def test_cm_column_matches_input(self, single_chr_map):
        tmp_dir, pattern = hapne_tmp_map(single_chr_map)
        try:
            df = read_shapeit(pattern.replace("@", "1"))
            assert list(df["cm"]) == pytest.approx([0.00, 0.10, 0.25, 0.50])
        finally:
            shutil.rmtree(tmp_dir)

    def test_rate_calculation(self, single_chr_map):
        """Rate = delta_cM / (delta_bp / 1e6), verified by hand.

        SNP1->SNP2: (0.10 - 0.00) / ((200000 - 100000) / 1e6) = 0.10 / 0.1 = 1.0 cM/Mb
        SNP2->SNP3: (0.25 - 0.10) / ((350000 - 200000) / 1e6) = 0.15 / 0.15 = 1.0 cM/Mb
        SNP3->SNP4: (0.50 - 0.25) / ((600000 - 350000) / 1e6) = 0.25 / 0.25 = 1.0 cM/Mb
        """
        tmp_dir, pattern = hapne_tmp_map(single_chr_map)
        try:
            df = read_shapeit(pattern.replace("@", "1"))
            assert list(df["rate"]) == pytest.approx([0.0, 1.0, 1.0, 1.0])
        finally:
            shutil.rmtree(tmp_dir)

    def test_rate_varies_correctly(self, tmp_path):
        """Verify rate reflects uneven spacing."""
        rows = [
            (1, "rs1", 0.0,  100_000),
            (1, "rs2", 1.0,  200_000),  # 1 cM over 0.1 Mb = 10 cM/Mb
            (1, "rs3", 1.5, 1_100_000), # 0.5 cM over 0.9 Mb = ~0.556 cM/Mb
        ]
        map_file = write_map(tmp_path, rows)
        tmp_dir, pattern = hapne_tmp_map(map_file)
        try:
            df = read_shapeit(pattern.replace("@", "1"))
            assert df["rate"].iloc[0] == pytest.approx(0.0)
            assert df["rate"].iloc[1] == pytest.approx(10.0)
            assert df["rate"].iloc[2] == pytest.approx(0.5 / 0.9, rel=1e-4)
        finally:
            shutil.rmtree(tmp_dir)


class TestHapneTmpMapMultiChrom:

    def test_one_file_per_chromosome(self, multi_chr_map):
        tmp_dir, pattern = hapne_tmp_map(multi_chr_map)
        try:
            for chrom in [1, 2, 3]:
                assert os.path.exists(pattern.replace("@", str(chrom)))
        finally:
            shutil.rmtree(tmp_dir)

    def test_no_extra_files(self, multi_chr_map):
        tmp_dir, pattern = hapne_tmp_map(multi_chr_map)
        try:
            files = os.listdir(tmp_dir)
            assert set(files) == {"chr1.txt", "chr2.txt", "chr3.txt"}
        finally:
            shutil.rmtree(tmp_dir)

    def test_each_chromosome_has_correct_row_count(self, multi_chr_map):
        tmp_dir, pattern = hapne_tmp_map(multi_chr_map)
        try:
            assert len(read_shapeit(pattern.replace("@", "1"))) == 3
            assert len(read_shapeit(pattern.replace("@", "2"))) == 2
            assert len(read_shapeit(pattern.replace("@", "3"))) == 2
        finally:
            shutil.rmtree(tmp_dir)

    def test_chromosomes_do_not_bleed_into_each_other(self, multi_chr_map):
        """Each file's cm values should match only that chromosome's input."""
        tmp_dir, pattern = hapne_tmp_map(multi_chr_map)
        try:
            df1 = read_shapeit(pattern.replace("@", "1"))
            df2 = read_shapeit(pattern.replace("@", "2"))
            assert list(df1["cm"]) == pytest.approx([0.00, 0.10, 0.30])
            assert list(df2["cm"]) == pytest.approx([0.00, 0.20])
        finally:
            shutil.rmtree(tmp_dir)

    def test_each_chromosome_first_rate_is_zero(self, multi_chr_map):
        tmp_dir, pattern = hapne_tmp_map(multi_chr_map)
        try:
            for chrom in [1, 2, 3]:
                df = read_shapeit(pattern.replace("@", str(chrom)))
                assert df["rate"].iloc[0] == 0.0, f"chr{chrom} first rate != 0"
        finally:
            shutil.rmtree(tmp_dir)


class TestHapneTmpMapEdgeCases:

    def test_duplicate_bp_no_exception(self, dup_bp_map):
        tmp_dir, pattern = hapne_tmp_map(dup_bp_map)
        shutil.rmtree(tmp_dir)  # just check it doesn't raise

    def test_duplicate_bp_rate_is_zero(self, dup_bp_map):
        tmp_dir, pattern = hapne_tmp_map(dup_bp_map)
        try:
            df = read_shapeit(pattern.replace("@", "1"))
            # The SNP at the duplicate position should have rate 0, not inf/nan
            assert df["rate"].iloc[1] == pytest.approx(0.0)
        finally:
            shutil.rmtree(tmp_dir)

    def test_duplicate_bp_no_nan_or_inf(self, dup_bp_map):
        tmp_dir, pattern = hapne_tmp_map(dup_bp_map)
        try:
            df = read_shapeit(pattern.replace("@", "1"))
            assert df["rate"].notna().all()
            assert np.isfinite(df["rate"]).all()
        finally:
            shutil.rmtree(tmp_dir)

    def test_single_snp_per_chromosome(self, single_snp_map):
        tmp_dir, pattern = hapne_tmp_map(single_snp_map)
        try:
            for chrom in [1, 2]:
                df = read_shapeit(pattern.replace("@", str(chrom)))
                assert len(df) == 1
                assert df["rate"].iloc[0] == 0.0
        finally:
            shutil.rmtree(tmp_dir)

    def test_cleanup_removes_tmp_dir(self, single_chr_map):
        tmp_dir, _ = hapne_tmp_map(single_chr_map)
        shutil.rmtree(tmp_dir)
        assert not os.path.exists(tmp_dir)

class TestHapneTmpMapRealData:

    MAP_FILE = os.path.join(os.path.dirname(__file__), "data", "genome_30chr.map")

    def test_produces_30_chromosome_files(self):
        tmp_dir, pattern = hapne_tmp_map(self.MAP_FILE)
        try:
            for chrom in range(1, 31):
                assert os.path.exists(pattern.replace("@", str(chrom)))
        finally:
            shutil.rmtree(tmp_dir)

    def test_no_extra_files(self):
        tmp_dir, pattern = hapne_tmp_map(self.MAP_FILE)
        try:
            files = os.listdir(tmp_dir)
            assert set(files) == {f"chr{c}.txt" for c in range(1, 31)}
        finally:
            shutil.rmtree(tmp_dir)

    def test_first_rate_zero_all_chromosomes(self):
        tmp_dir, pattern = hapne_tmp_map(self.MAP_FILE)
        try:
            for chrom in range(1, 31):
                df = read_shapeit(pattern.replace("@", str(chrom)))
                assert df["rate"].iloc[0] == 0.0, f"chr{chrom} first rate != 0"
        finally:
            shutil.rmtree(tmp_dir)

    def test_rates_approximately_one_cM_per_Mb(self):
        """All non-first SNPs should have rate ~1.0 cM/Mb given 1cM = 1Mb."""
        tmp_dir, pattern = hapne_tmp_map(self.MAP_FILE)
        try:
            for chrom in range(1, 31):
                df = read_shapeit(pattern.replace("@", str(chrom)))
                non_first_rates = df["rate"].iloc[1:]
                assert non_first_rates.between(0.5, 2.0).all(), \
                    f"chr{chrom} has rates outside expected range: {non_first_rates.describe()}"
        finally:
            shutil.rmtree(tmp_dir)

    def test_no_nan_or_inf_in_any_chromosome(self):
        tmp_dir, pattern = hapne_tmp_map(self.MAP_FILE)
        try:
            for chrom in range(1, 31):
                df = read_shapeit(pattern.replace("@", str(chrom)))
                assert df["rate"].notna().all(), f"chr{chrom} has NaN rates"
                assert np.isfinite(df["rate"]).all(), f"chr{chrom} has inf rates"
        finally:
            shutil.rmtree(tmp_dir)

    def test_cm_values_preserved(self):
        """cm column in output should match input cm values exactly."""
        input_df = pd.read_csv(
            self.MAP_FILE, sep="\t", header=None, names=["chrom", "snp_id", "cm", "bp"]
        )
        tmp_dir, pattern = hapne_tmp_map(self.MAP_FILE)
        try:
            for chrom in range(1, 31):
                out_df = read_shapeit(pattern.replace("@", str(chrom)))
                expected_cm = input_df[input_df["chrom"] == chrom]["cm"].values
                assert out_df["cm"].values == pytest.approx(expected_cm)
        finally:
            shutil.rmtree(tmp_dir)

    def test_bp_values_preserved(self):
        """bp column in output should match input bp values exactly."""
        input_df = pd.read_csv(
            self.MAP_FILE, sep="\t", header=None, names=["chrom", "snp_id", "cm", "bp"]
        )
        tmp_dir, pattern = hapne_tmp_map(self.MAP_FILE)
        try:
            for chrom in range(1, 31):
                out_df = read_shapeit(pattern.replace("@", str(chrom)))
                expected_bp = input_df[input_df["chrom"] == chrom]["bp"].values
                assert list(out_df["bp"]) == list(expected_bp)
        finally:
            shutil.rmtree(tmp_dir)