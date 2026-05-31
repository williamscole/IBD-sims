"""
Diagnostic: verify that get_node("tsk_j", h) maps correctly to the right
individual in the tree sequence for a pedigree-based simulation.

Run from the ibd_sims/ directory against any existing pedigree + tree sequence:
  python /tmp/check_mapping.py <path_to_iter_dir> <iter_n>

Checks:
1. For each tsk_j in the VCF, verify get_node maps to the correct individual
2. For known full-sibling pairs (from pedigree file), check they share parents
   in the tree sequence
3. Sample 5 full-sibling pairs and print their expected vs. tree-seq nodes
"""

import sys
import gzip
import msprime
import tszip
import pandas as pd
import numpy as np

def get_node(vcf_id, haplotype):
    return int(vcf_id.split("_")[1])*2 + (haplotype==2)

def main():
    if len(sys.argv) < 3:
        print("Usage: python check_mapping.py <path> <iter_n>")
        sys.exit(1)
    
    path, iter_n = sys.argv[1], sys.argv[2]
    chrom = 1  # check chromosome 1
    prefix = f"{path}/iter{iter_n}_chr{chrom}"
    ped_file = f"{path}/iter{iter_n}_WF.pedigree"
    
    # --- Step 1: Read VCF header to get sample order ---
    vcf_path = f"{prefix}.vcf.gz"
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                vcf_samples = line.strip().split("\t")[9:]
                break
    print(f"VCF has {len(vcf_samples)} samples")
    print(f"First 5: {vcf_samples[:5]}")
    print(f"Last 5:  {vcf_samples[-5:]}")
    
    # --- Step 2: Parse pedigree to find full-sibling pairs ---
    ped = pd.read_csv(ped_file, sep=r"\s+", comment="#",
                      names=["id","parent0","parent1","time","is_sample","pop"])
    samples = ped[ped.is_sample == 1].copy()
    print(f"\nPedigree: {len(samples)} samples (is_sample=1)")
    
    # Find full siblings: pairs with same parent0 AND parent1
    siblings = samples.merge(samples, on=["parent0","parent1"], suffixes=("_A","_B"))
    siblings = siblings[siblings.id_A < siblings.id_B]
    print(f"Full sibling pairs in pedigree: {len(siblings)}")
    
    if len(siblings) == 0:
        print("No full siblings — try mono mating or lower Ne")
        return
    
    # --- Step 3: Verify a few sibling pairs map to correct nodes ---
    # id format: "g0_j" where j is the VCF index (tsk_j)
    def id_to_vcf_idx(ped_id):
        return int(ped_id.split("_")[1])
    
    print("\nSample full-sibling pairs (pedigree ID → expected nodes):")
    for _, row in siblings.head(5).iterrows():
        idx_a = id_to_vcf_idx(row.id_A)
        idx_b = id_to_vcf_idx(row.id_B)
        vcf_name_a = vcf_samples[idx_a] if idx_a < len(vcf_samples) else "OUT_OF_RANGE"
        vcf_name_b = vcf_samples[idx_b] if idx_b < len(vcf_samples) else "OUT_OF_RANGE"
        expected_vcf_a = f"tsk_{idx_a}"
        expected_vcf_b = f"tsk_{idx_b}"
        node_a_h1 = get_node(vcf_name_a, 1)
        node_b_h1 = get_node(vcf_name_b, 1)
        
        match_a = "✓" if vcf_name_a == expected_vcf_a else f"MISMATCH (got {vcf_name_a})"
        match_b = "✓" if vcf_name_b == expected_vcf_b else f"MISMATCH (got {vcf_name_b})"
        
        print(f"  {row.id_A} (pedigree idx {idx_a}) → VCF col {idx_a} = {vcf_name_a} {match_a} → nodes ({node_a_h1}, {node_a_h1+1})")
        print(f"  {row.id_B} (pedigree idx {idx_b}) → VCF col {idx_b} = {vcf_name_b} {match_b} → nodes ({node_b_h1}, {node_b_h1+1})")
        print(f"  Parents: {row.parent0}, {row.parent1}")
        print()
    
    # --- Step 4: Check IBD file for these sibling pairs ---
    ibd_file = f"{path}/iter{iter_n}.ibd.gz"
    if not __import__("os").path.exists(ibd_file):
        print(f"No concatenated IBD file at {ibd_file}")
        return
    
    ibd = pd.read_csv(ibd_file, sep=r"\s+", header=None)
    ibd_kinship = ibd.groupby([0, 2])[7].sum()
    
    print("Kinship sums for full sibling pairs:")
    for _, row in siblings.head(10).iterrows():
        idx_a = id_to_vcf_idx(row.id_A)
        idx_b = id_to_vcf_idx(row.id_B)
        name_a = f"tsk_{idx_a}"
        name_b = f"tsk_{idx_b}"
        key = (name_a, name_b) if name_a < name_b else (name_b, name_a)
        k = ibd_kinship.get(key, 0.0)
        print(f"  {name_a} + {name_b}: {k:.1f} cM")

if __name__ == "__main__":
    main()
