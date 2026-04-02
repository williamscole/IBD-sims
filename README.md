# ibdne-sims

A pipeline for simulating realistic IBD segments under a range of demographic histories and mating models, designed for evaluating [IBDNe](https://faculty.washington.edu/browning/ibdne.html) (a tool that infers effective population size from IBD segments). The simulated IBD segments can also be used for any other downstream analysis.

Ancestry is simulated with [msprime](https://tskit.dev/msprime/docs/stable/intro.html), supporting both standard coalescent simulations and simulations that begin with an explicit Wright-Fisher pedigree for recent generations (the DTWF model). IBD segments are detected from simulated genotype data using [hap-ibd](https://github.com/browning-lab/hap-ibd).

## Requirements

**Python packages:** msprime, numpy, pandas, matplotlib, seaborn, PyYAML, submitit, stdpopsim

**External tools (must be on `$PATH` or configured in `setup.yaml`):**

- [hap-ibd](https://github.com/browning-lab/hap-ibd) — IBD segment detection (Java jar)
- [IBDNe](https://faculty.washington.edu/browning/ibdne.html) — Ne inference from IBD (Java jar)
- Java (JRE 8+)
- bcftools (only required if `keep_all_files: true`)

## Setup

Before running anything, edit `setup.yaml` to point to your local paths:

```yaml
maf_pickle: ukb_snps.pkl                    # SNP density/MAF pickle (included in repo)
hap_ibd_jar: /path/to/hap-ibd.jar           # hap-ibd Java jar
ibdne_jar: /path/to/ibdne.jar               # IBDNe Java jar
genetic_map_chr1: /path/to/genetic_map_GRCh37_chr1.txt.gz  # HapMap genetic map (only needed for 22-chromosome sims)
```

## Quick start

The main entry point is `simulate.py`. It handles pedigree creation, per-chromosome simulation, post-processing (concatenation, purple node matrices, TMRCA tracking), IBDNe inference, and plotting — all in one command.

**Run on a Slurm cluster (default):**

```bash
python simulate.py yaml_files/arg1.yaml
```

**Run locally (no Slurm):**

```bash
python simulate.py yaml_files/arg1.yaml --local --workers 8
```

**Test with the debug config first** to make sure everything is wired up:

```bash
python simulate.py yaml_files/debug.yaml --local
```

This runs a single replicate over 2 small chromosomes using a pure coalescent model — it should finish in a few minutes.

**Resume a previous run** by passing the output directory instead of a YAML file:

```bash
python simulate.py path/to/existing/run/
```

The pipeline detects which jobs already completed and only re-runs what's missing.

**Override YAML parameters on the command line** without editing the file:

```bash
python simulate.py yaml_files/arg1.yaml --set iter=5 --set pedigree.mating=mono --set run_ibdne=false
```

## What `simulate.py` does

The pipeline runs in three phases:

1. **Pedigree creation** (if `pedigree.pedigree_mode: true`): generates a Wright-Fisher pedigree file for each replicate using `wf_pedigree.py`.

2. **Simulation** (parallelized across replicates × chromosomes): for each (replicate, chromosome) pair, `simulations.py` runs msprime ancestry simulation, adds mutations, subsets to realistic SNP density/MAF, writes a VCF + genetic map, runs hap-ibd to detect IBD segments, and computes TMRCAs for each segment.

3. **Post-processing** (per replicate): concatenates per-chromosome IBD, map, and TMRCA files; computes purple node matrices; optionally runs IBDNe; and generates a summary plot comparing inferred Ne to the true demographic history.

## Output structure

After a run completes, the output directory contains:

```
output_dir/
├── args.yaml              # Copy of the configuration used
├── iter1_chr1.ibd.gz      # Per-chromosome IBD (deleted after concatenation)
├── iter1.ibd.gz           # Concatenated IBD segments for replicate 1
├── iter1.map              # Concatenated genetic map
├── iter1.tmrca.gz         # TMRCA annotations for each IBD segment
├── iter1.npy              # Purple node matrix
├── ibdne/                 # IBDNe results (or dir_name from config)
│   ├── args.yaml
│   ├── iter1.ne           # IBDNe output
│   └── ...
├── plot.png               # Ne estimates vs truth
├── slurm/                 # Slurm/submitit logs
└── errors/                # Error logs for failed jobs
```

## Simulation configuration

Each simulation is configured via a YAML file. A template with default values and comments is provided in `yaml_files/args.yaml`. The key parameters are:

**Simulation**

- `iter` — number of independent replicates
- `samples` — number of diploid individuals
- `end_chr` — number of chromosomes. `30` simulates 30 chromosomes of 100 Mb each with a flat recombination rate. `22` simulates the 22 human autosomes using the HapMap GRCh37 genetic map (path set in `setup.yaml`)
- `base_dir` — optional parent directory for output (default: current working directory)
- `gb` — memory allocation in GB (used for Java heap and Slurm `--mem`)
- `sim_time` — Slurm wall-time for simulation jobs (e.g., `--time=1:30:00`)

**Demographic model**

- `custom_demo.path` — Python file containing an `msprime.Demography` object
- `custom_demo.object` — name of the variable in that file

Several models are provided in `demography.py`: `constant_Ne` (10k), `constant_Ne100k` (100k), `euro_bottleneck`, `himba`, `ooa2` (two-population Out-of-Africa), and `ashkenazi` (via stdpopsim).

**Mating model**

- `pedigree.pedigree_mode` — if `true`, recent generations use an explicit Wright-Fisher pedigree before switching to the coalescent
- `pedigree.mating` — `di` (random mating) or `mono` (monogamous)
- `pedigree.gen_end` — number of pedigree generations before handing off to the coalescent

**IBDNe**

- `run_ibdne` — whether to run IBDNe on the detected IBD segments
- `ibdne_time` — Slurm wall-time for IBDNe jobs
- `nboots`, `nits`, `mincm`, `trimcm`, `gmin`, `gmax` — IBDNe parameters (see IBDNe documentation)
- `ibd_filter` — optional IBD filtering strategy before running IBDNe: `null` (no filtering), `related` (subset to close relatives), or `unrelated` (prune relatives)

## Example configurations

Ready-to-run configurations are provided in `yaml_files/`:

| File | Demographic model | Samples | Mating | Pedigree gens | Replicates |
|------|-------------------|---------|--------|---------------|------------|
| `debug.yaml` | Constant Ne (10k) | 1,000 | coalescent | — | 1 |
| `arg1.yaml` | Constant Ne (10k) | 1,000 | Random | 25 | 50 |
| `arg2.yaml` | Constant Ne (10k) | 1,000 | Monogamous | 25 | 50 |
| `arg3.yaml` | Constant Ne (100k) | 1,000 | Random | 25 | 50 |
| `arg4.yaml` | Constant Ne (100k) | 1,000 | Monogamous | 25 | 50 |
| `arg5.yaml` | Out-of-Africa (2-pop) | 2,000 | Random | 25 | 50 |
| `arg6.yaml` | Out-of-Africa (2-pop) | 2,000 | Monogamous | 25 | 50 |
| `arg7.yaml` | Quebec (empirical) | 10,000 | — | — | 1 |
| `arg8.yaml` | Ashkenazi | 1,000 | Random | 12 | 50 |

## SNP density and MAF distribution

Simulated VCFs are thinned to match a realistic SNP density and minor allele frequency distribution. This is controlled by a pickle file (`ukb_snps.pkl`) that stores MAF bin proportions and mean SNP density, used in `write_vcf.py` to subset mutations from the simulated tree sequence.

A default `ukb_snps.pkl` derived from UK Biobank genotype data is included. To generate your own from a different dataset:

```bash
python maf_buckets.py \
    --afreq-chr1 /path/to/chr1.afreq \
    --bim-chr1 /path/to/chr1.bim \
    --output my_snps.pkl
```

`--afreq-chr1` should be `plink2 --freq` output for chromosome 1; other chromosomes are found by replacing `chr1` with `chr{n}`. Same convention for `--bim-chr1`. Then set `maf_pickle` in `setup.yaml` to point to your output file.

## Analysis tools

`analysis_tools.py` provides post-hoc analysis functions that can be run on completed simulation outputs:

```bash
# Compute TMRCA quantiles by IBD segment length
python analysis_tools.py tmrca_age path/to/run/

# Compute purple node matrices from IBD files
python analysis_tools.py purple_matrix path/to/run/

# Compute TMRCA quantiles stratified by kinship degree
python analysis_tools.py tmrca_kinship path/to/run/
```

## Adding a new demographic model

1. Define an `msprime.Demography` object in a Python file (e.g., `demography.py` or your own file).
2. Reference it in a YAML config:

```yaml
custom_demo:
  path: demography.py
  object: my_new_model
```

3. Run `python simulate.py yaml_files/my_config.yaml`.

## File overview

| File | Purpose |
|------|---------|
| `simulate.py` | Main orchestrator — run this |
| `simulations.py` | Core simulation logic (msprime, VCF, hap-ibd, TMRCA) |
| `demography.py` | Demographic model definitions |
| `wf_pedigree.py` | Wright-Fisher pedigree generation |
| `write_vcf.py` | VCF and genetic map output with realistic SNP thinning |
| `maf_buckets.py` | Build SNP density/MAF pickle from plink files |
| `plot.py` | Plot IBDNe estimates vs true Ne |
| `filter_ibd.py` | IBD filtering (related/unrelated subsets) |
| `purple.py` | Purple node matrix computation |
| `concat_tmrca.py` | Concatenate per-chromosome TMRCA files |
| `kinship.py` | Kinship degree classification from IBD |
| `analysis_tools.py` | Post-hoc analysis (TMRCA by length, kinship, etc.) |
| `yaml_tools.py` | Utilities for generating/modifying YAML configs in batch |
| `setup.yaml` | Path configuration — edit before first run |
| `yaml_files/args.yaml` | Default YAML template with all parameters |
| `yaml_files/arg*.yaml` | Pre-built configurations |