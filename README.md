# ibdne-sims

A pipeline for generating realistic IBD segments under a range of demographic histories and mating models. Ancestry is simulated using [msprime](https://tskit.dev/msprime/docs/stable/intro.html), supporting both standard coalescent simulations and simulations that begin with an explicit Wright-Fisher pedigree for recent generations (DTWF model). IBD segments are detected from simulated genotype data using [hap-ibd](https://github.com/browning-lab/hap-ibd).

The pipeline was developed to evaluate [IBDNe](https://faculty.washington.edu/browning/ibdne.html), a tool that infers effective population size (Ne) through time from IBD segments, but the simulated IBD segments can be used for any downstream purpose.

## SNP density and MAF distribution

Simulated VCFs are generated to match a realistic SNP density and minor allele frequency (MAF) distribution. This is controlled by a pickle file (`ukb_snps.pkl`) that stores the MAF bin proportions and mean SNP density used to subset mutations from the simulated tree sequence.

A default `ukb_snps.pkl` derived from UK Biobank genotype data is included in the repository. If you would like to generate your own from a different dataset, you can do so using `maf_buckets.py`:

```bash
python maf_buckets.py \
    --afreq-chr1 /path/to/chr1.afreq \
    --bim-chr1 /path/to/chr1.bim \
    --output ukb_snps.pkl
```

`--afreq-chr1` should be the output of `plink2 --freq` for chromosome 1. Files for other chromosomes are expected in the same directory with `chr1` replaced by `chr{n}`. The same convention applies to `--bim-chr1`. The output path passed to `--output` should then be set as `maf_pickle` in `setup.yaml`.

## Simulation configuration

Each simulation is configured via a YAML file. A template with default values is provided in `yaml_files/args.yaml`. The key parameters are:

**Simulation**
- `iter` — number of independent simulation replicates to run
- `samples` — number of diploid individuals to simulate
- `end_chr` — number of chromosomes to simulate. If `30`, simulates 30 chromosomes of 100 Mb each with a flat recombination rate of 1e-8. If `22`, simulates the 22 autosomes using the real HapMap GRCh37 genetic map (path set in `setup.yaml`)
- `base_dir` — optional base directory under which the output directory will be created. If not set, the output directory is created in the current working directory

**Demographic model**
- `custom_demo` — specifies the demographic history to simulate under. Set `path` to a Python file and `object` to the name of an `msprime.Demography` object defined in that file. Several models are provided in `demography.py`

**Mating model**
- `pedigree.pedigree_mode` — if `true`, recent generations are simulated using an explicit Wright-Fisher pedigree before switching to the coalescent
- `pedigree.mating` — mating model for the pedigree phase; either `di` (random mating) or `mono` (monogamous)
- `pedigree.gen_end` — number of generations to simulate under the pedigree model before handing off to the coalescent

**IBDNe**
- `run_ibdne` — whether to run IBDNe on the simulated IBD segments

## Example configurations

Eight ready-to-run configurations are provided in `yaml_files/`:

| File | Demographic model | N samples | Mating | Pedigree generations |
|------|------------------|-----------|--------|----------------------|
| `debug.yaml` | Constant Ne (10k) | 1,000 | — | — |
| `arg1.yaml` | Constant Ne (10k) | 1,000 | Random | 25 |
| `arg2.yaml` | Constant Ne (10k) | 1,000 | Monogamous | 25 |
| `arg3.yaml` | Constant Ne (100k) | 1,000 | Random | 25 |
| `arg4.yaml` | Constant Ne (100k) | 1,000 | Monogamous | 25 |
| `arg5.yaml` | Out-of-Africa (2-pop) | 2,000 | Random | 25 |
| `arg6.yaml` | Out-of-Africa (2-pop) | 2,000 | Monogamous | 25 |
| `arg7.yaml` | Quebec (real data) | 10,000 | — | — |
| `arg8.yaml` | Ashkenazi | 1,000 | Random | 12 |

All configurations except `arg_debug.yaml` use 50 replicates and 30 simulated chromosomes of 100 Mb each. The debug configuration runs a single replicate over 2 chromosomes using a pure coalescent simulation and is intended for testing the pipeline end-to-end. The Quebec configuration uses empirical tree sequences rather than a simulated demographic model.