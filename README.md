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
hapmap_chr1: /path/to/genetic_map_GRCh37_chr1.txt.gz
```

`hapmap_chr1` is only needed when `end_chr: 22` (real autosomes). For simulated chromosomes (`end_chr: 30`) it is ignored.

## Quick start

The main entry point for simulation is `simulate.py`. It handles pedigree creation, per-chromosome simulation, and concatenation of output files.

**Test with the debug config first** to make sure everything is wired up:

```bash
python simulate.py yaml_files/debug.yaml --local
```

**Run on a Slurm cluster (default):**

```bash
python simulate.py yaml_files/arg1.yaml
```

**Run locally:**

```bash
python simulate.py yaml_files/arg1.yaml --local --workers 8
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

The pipeline has two independent phases.

### Phase 1: Simulation

`simulate.py` produces IBD segments. For each replicate, the outputs are:

- `iter{n}.ibd.gz` — concatenated IBD segments across all chromosomes
- `iter{n}.tmrca.gz` — TMRCA annotations for each IBD segment
- `iter{n}.map` — concatenated genetic map

Optional outputs (with `keep_all_files: true`): concatenated VCF, HBD, and per-chromosome tree sequences.

### Phase 2: Post-processing

`post_process.py` runs analyses on the simulation output. It can be run independently of the simulation, and re-run with different parameters without re-simulating.

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

Each simulation is configured via a YAML file with two sections: simulation parameters at the top level, and post-processing module blocks below.

### Simulation parameters

```yaml
# Housekeeping
base_dir: ibd_sims          # parent directory for output (null = current dir)
dir_name: ibd-sims          # subdirectory name within the run directory
keep_all_files: false       # if true, also saves VCF and HBD files

# Computational resources (simulation)
gb: 8                       # memory in GB for simulation jobs
sim_min: 30                 # wall-time in minutes for simulation jobs
nthreads: 8                 # threads for hap-ibd

# Simulation
iter: 50                    # number of independent replicates
samples: 1000               # number of diploid individuals
end_chr: 30                 # chromosomes to simulate (30 = 30x100Mb; 22 = real autosomes)

# Demographic model (one of custom_demo or custom_sim must be set)
custom_demo:
  path: demography.py       # Python file with an msprime.Demography object
  object: constant_Ne       # variable name in that file
custom_sim:
  path: null                # Python file with a custom tree-sequence loader
  object: null

# Mating model
pedigree:
  pedigree_mode: true       # use Wright-Fisher pedigree for recent generations
  mating: di                # "di" (random) or "mono" (monogamous)
  gen_end: 25               # pedigree generations before coalescent takeover
  pedigree_file: null       # path to pre-existing pedigree file (optional)
```

### Post-processing parameters

Set `post_process` to a comma-separated list of modules to run, or `null` to skip all post-processing. Each module has its own nested config block with `object` (Python class name), `path` (Python file), analysis parameters, and optional resource overrides.

Resource fields (`workers`, `mem_gb`, `time_min`) can be set at the top level as defaults and overridden per module:

```yaml
post_process: ibdne,purple_nodes  # null to skip; options: ibdne,hapne_ibd,hapne_ld,purple_nodes

# Default resources for all post-processing jobs
workers: 8
mem_gb: 16
time_min: 120

ibdne:
  object: PostProcessIBDNe
  path: post_modules.py
  ibd_filter: null
  filtersamples: false
  mincm: 2
  trimcm: 0.2
  gmin: 1
  gmax: 300
  nboots: 80
  nits: 1000
  npairs: 0
  workers: 8       # module-level override
  mem_gb: 16
  time_min: 120

hapne_ibd:
  object: PostProcessHapNeIBD
  path: post_modules.py
  filter: unfiltered
  workers: 4
  mem_gb: 16
  time_min: 120

hapne_ld:
  object: PostProcessHapNeLD
  path: post_modules.py
  filter: unfiltered
  workers: 4
  mem_gb: 16
  time_min: 120

purple_nodes:
  object: PostProcessPurple
  path: post_modules.py
  workers: 4
  mem_gb: 8
  time_min: 60
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

All configs have `post_process: null` by default. To run post-processing, edit the file or use `--set`:

```bash
python simulate.py yaml_files/arg1.yaml --set post_process=ibdne
```

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

## Adding a new demographic model

Define an `msprime.Demography` object in a Python file and reference it in the YAML:

```yaml
custom_demo:
  path: my_demography.py
  object: my_model
```

Several models are provided in `demography.py`: `constant_Ne` (10k), `constant_Ne100k` (100k), `euro_bottleneck`, `himba`, `ooa2` (two-population Out-of-Africa), and `ashkenazi` (via stdpopsim).

## Writing a custom post-processing module

Subclass `PostProcessor` from `post_process.py`:

```python
# my_analysis.py
from post_process import PostProcessor

class MyAnalysis(PostProcessor):
    sub_config_key = "my_analysis"
    resource_fields = ["local", "workers", "mem_gb", "time_min"]

    def execute(self):
        self._execute_helper()
        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop()

    def _single_iter(self, iter_n):
        cfg = self._get_sub_config()
        prefix = f"{self.path}/iter{iter_n}"
        # Read from: {prefix}.ibd.gz, {prefix}.map, {prefix}.tmrca.gz
        # Write to:  self.out_dir
```

Then add it to your YAML:

```yaml
post_process: ibdne,my_analysis

my_analysis:
  object: MyAnalysis
  path: my_analysis.py
  my_param: 42
  workers: 4
  time_min: 60
```

Key things to know:

- `sub_config_key` must match the YAML block name exactly.
- `resource_fields` lists fields excluded from config comparison when deciding whether to reuse an existing output directory.
- `self._get_sub_config()` returns the config object for your module, with all YAML values as attributes.
- `self._get_resource(name)` checks the module config first, then falls back to top-level defaults. Use this for `workers`, `mem_gb`, and `time_min`.
- `self._execute_loop()` handles both local and Slurm execution based on the `local` resource setting.
- `self.out_dir` is set by `_execute_helper()` and points to the numbered output subdirectory.
- Check `self.single_iter` to distinguish single-iteration vs. full runs.

## SNP density and MAF distribution

Simulated VCFs are thinned to match a realistic SNP density and minor allele frequency distribution, controlled by `ukb_snps.pkl`. A default pickle derived from UK Biobank genotype data is included. To generate your own:

```bash
python maf_buckets.py \
    --afreq-chr1 /path/to/chr1.afreq \
    --bim-chr1 /path/to/chr1.bim \
    --output my_snps.pkl
```

Then point `maf_pickle` in `setup.yaml` to your output file.

## File overview

| File | Purpose |
|------|---------|
| `simulate.py` | Simulation orchestrator — entry point for Phase 1 |
| `simulations.py` | Core simulation logic (msprime, VCF, hap-ibd, TMRCA) |
| `post_process.py` | Post-processing orchestrator — entry point for Phase 2, `PostProcessor` ABC |
| `post_modules.py` | Built-in post-processors: IBDNe, HapNe-IBD, HapNe-LD, purple nodes |
| `demography.py` | Demographic model definitions |
| `wf_pedigree.py` | Wright-Fisher pedigree generation |
| `write_vcf.py` | VCF and genetic map output with realistic SNP thinning |
| `maf_buckets.py` | Build SNP density/MAF pickle from plink files |
| `run_hapne.py` | HapNe-IBD and HapNe-LD runner utilities |
| `plot.py` | Plot IBDNe estimates vs true Ne |
| `filter_ibd.py` | IBD filtering (related/unrelated subsets) |
| `purple.py` | Purple node matrix computation |
| `concat_tmrca.py` | Concatenate per-chromosome TMRCA files |
| `kinship.py` | Kinship degree classification from IBD |
| `analysis_tools.py` | Post-hoc analysis (TMRCA by length, kinship, etc.) |
| `setup.yaml` | Machine-specific path configuration — edit before first run |
| `yaml_files/` | Per-experiment simulation configs |