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
maf_pickle: ukb_snps.pkl
hap_ibd_jar: /path/to/hap-ibd.jar
ibdne_jar: /path/to/ibdne.jar
genetic_map_chr1: /path/to/genetic_map_GRCh37_chr1.txt.gz
```

## Quick start

The main entry point for simulation is `simulate.py`. It handles pedigree creation, per-chromosome simulation, and concatenation of output files.

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

**Resume a previous run** by passing the output directory instead of a YAML file:

```bash
python simulate.py path/to/existing/run/
```

**Override YAML parameters on the command line:**

```bash
python simulate.py yaml_files/arg1.yaml --set iter=5 --set pedigree.mating=mono
```

## Pipeline overview

The pipeline has two independent phases:

### Phase 1: Simulation

`simulate.py` produces IBD segments. For each replicate, the outputs are:

- `iter{n}.ibd.gz` — concatenated IBD segments across all chromosomes
- `iter{n}.tmrca.gz` — TMRCA annotations for each IBD segment
- `iter{n}.map` — concatenated genetic map

Optional outputs (with `keep_all_files: true`): concatenated VCF, HBD, and per-chromosome tree sequences.

### Phase 2: Post-processing

`post_process.py` runs analyses on the simulation output. It can be run independently of the simulation, and can be re-run with different parameters without re-simulating.

```bash
# Run post-processing on a completed simulation
python post_process.py path/to/run/

# Re-run IBDNe with different parameters
python post_process.py path/to/run/ --set ibdne.nboots=100 ibdne.mincm=3

# Run via Slurm instead of locally
python post_process.py path/to/run/ --set local=false
```

Post-processing output goes into numbered subdirectories (e.g., `ibdne/001/`, `ibdne/002/`). If you re-run with the same analysis parameters, it overwrites the matching directory. If you change parameters, it creates a new one. Each subdirectory contains an `args.yaml` recording exactly which parameters were used for that run.

## YAML configuration

Each simulation is configured via a YAML file. The file has two sections: simulation parameters at the top level, and post-processing modules as nested blocks.

### Simulation parameters

```yaml
# Simulation
iter: 50                    # number of independent replicates
samples: 1000               # number of diploid individuals
end_chr: 30                 # chromosomes to simulate (30 = 30x100Mb; 22 = real autosomes)
sim_time: --time=1:30:00    # Slurm wall-time for simulation jobs
base_dir: null              # optional parent directory for output

# Demographic model
custom_demo:
  path: demography.py       # Python file with msprime.Demography object
  object: constant_Ne       # variable name in that file

# Mating model
pedigree:
  pedigree_mode: true       # use Wright-Fisher pedigree for recent generations
  mating: di                # "di" (random) or "mono" (monogamous)
  gen_end: 25               # pedigree generations before coalescent takeover
```

### Post-processing parameters

The `post_process` field is a comma-separated list of analysis modules to run (or `none`). Each module has its own nested configuration block, which must include `object` (the Python class name) and `path` (the Python file containing it).

```yaml
post_process: ibdne,purple_nodes   # or "none"

ibdne:
  object: PostProcessIBDNe
  path: post_modules.py
  ibd_filter: null
  mincm: 2
  nboots: 80
  nits: 1000
  npairs: 0
  filtersamples: false
  trimcm: 0.2
  gmin: 1
  gmax: 300
  workers: 8                # resource override (top-level default: 1)

purple_nodes:
  object: PostProcessPurple
  path: post_modules.py
```

Resource fields (`local`, `workers`, `mem_gb`, `time_min`) can be set at the top level as defaults, or inside any module block to override for that specific analysis. Module-level values take precedence.

```yaml
# Top-level defaults
local: true
workers: 1
mem_gb: 8
time_min: 60

ibdne:
  workers: 8      # IBDNe gets 8 CPUs; other modules use the default of 1
  mem_gb: 16      # IBDNe gets 16 GB; others get 8
```

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

## Writing a custom post-processing module

You can add your own analysis step by subclassing `PostProcessor` from `post_process.py`. Here's the minimal pattern:

```python
# my_analysis.py
from post_process import PostProcessor

class MyAnalysis(PostProcessor):
    sub_config_key = "my_analysis"        # matches the YAML block name
    resource_fields = ["local", "workers", "mem_gb", "time_min"]  # excluded from config comparison

    def execute(self):
        self._execute_helper()            # creates output dir, dumps config

        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop()          # runs _single_iter across all iterations

    def _single_iter(self, iter_n):
        cfg = self._get_sub_config()      # access my_analysis.* config values
        prefix = f"{self.path}/iter{iter_n}"

        # Your analysis logic here.
        # Read from: {prefix}.ibd.gz, {prefix}.map, {prefix}.tmrca.gz
        # Write to:  self.out_dir
```

Then add it to your YAML:

```yaml
post_process: ibdne,my_analysis

my_analysis:
  object: MyAnalysis
  path: my_analysis.py
  my_param: 42               # accessible as cfg.my_param
  workers: 4                  # resource override
```

Key things to know:

- `sub_config_key` must match the YAML block name exactly.
- `resource_fields` lists attribute names that should be ignored when comparing configs to decide whether to reuse an existing output directory. Typically these are infrastructure settings that don't affect results.
- `self._get_sub_config()` returns the `PostProcessConfig` object for your module, with all YAML values set as attributes.
- `self._get_resource(name)` checks the module config first, then falls back to top-level defaults. Use this for `workers`, `mem_gb`, `local`, and `time_min`.
- `self._execute_loop()` handles both local (sequential or `ProcessPoolExecutor`) and Slurm (`submitit`) execution based on the `local` resource setting.
- `self.out_dir` is set by `_execute_helper()` and points to the numbered output subdirectory (e.g., `path/my_analysis/001/`).
- The processor can run on all iterations (`n_iter` is set) or a single iteration (`iter_n` is set), but not both. Check `self.single_iter` to distinguish.

## SNP density and MAF distribution

Simulated VCFs are thinned to match a realistic SNP density and minor allele frequency distribution. This is controlled by `ukb_snps.pkl`, which stores MAF bin proportions and mean SNP density.

A default pickle derived from UK Biobank genotype data is included. To generate your own:

```bash
python maf_buckets.py \
    --afreq-chr1 /path/to/chr1.afreq \
    --bim-chr1 /path/to/chr1.bim \
    --output my_snps.pkl
```

Then set `maf_pickle` in `setup.yaml` to point to your output file.

## Adding a new demographic model

Define an `msprime.Demography` object in a Python file and reference it in the YAML:

```yaml
custom_demo:
  path: my_demography.py
  object: my_model
```

Several models are provided in `demography.py`: `constant_Ne` (10k), `constant_Ne100k` (100k), `euro_bottleneck`, `himba`, `ooa2` (two-population Out-of-Africa), and `ashkenazi` (via stdpopsim).

## Output structure

```
run_directory/
├── args.yaml                 # simulation config
├── iter1.ibd.gz              # concatenated IBD segments
├── iter1.map                 # concatenated genetic map
├── iter1.tmrca.gz            # TMRCA annotations
├── ibdne/
│   ├── 001/
│   │   ├── args.yaml         # analysis config for this run
│   │   ├── iter1.ne          # IBDNe output
│   │   └── ...
│   └── 002/                  # re-run with different parameters
│       ├── args.yaml
│       └── ...
├── purple_nodes/
│   └── 001/
│       ├── args.yaml
│       ├── iter1.npy
│       └── ...
├── slurm/                    # Slurm/submitit logs
└── errors/                   # error logs
```

## File overview

| File | Purpose |
|------|---------|
| `simulate.py` | Phase 1 orchestrator — simulation entry point |
| `simulations.py` | Core simulation logic (msprime, VCF, hap-ibd, TMRCA) |
| `post_process.py` | Phase 2 orchestrator — post-processing entry point, `PostProcessor` ABC |
| `post_modules.py` | Built-in post-processors: IBDNe, purple nodes |
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
| `setup.yaml` | Path configuration — edit before first run |
| `yaml_files/args.yaml` | Default YAML template with all parameters |