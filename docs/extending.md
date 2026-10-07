# Extending the pipeline

- [Adding a demographic model](#adding-a-demographic-model)
- [Supplying your own tree sequences](#supplying-your-own-tree-sequences)
- [Writing a custom analysis](#writing-a-custom-analysis)
- [Using your own SNP density and allele frequencies](#using-your-own-snp-density-and-allele-frequencies)
- [How a simulation runs](#how-a-simulation-runs)
- [Repository layout](#repository-layout)

## Adding a demographic model

Define an [`msprime.Demography`](https://tskit.dev/msprime/docs/stable/demography.html) object in any Python file and point your config at it:

```yaml
custom_demo:
  path: my_demography.py
  object: my_model
```

Samples are drawn from the population named `pop_0`, and the "true Ne" in plots is that population's size.

Models already in `ibd_sims/demography.py`:

| Name | Description |
|------|-------------|
| `constant_Ne` | Constant Ne of 10,000 |
| `constant_Ne100k` | Constant Ne of 100,000 |
| `euro_bottleneck` | European-like: constant Ne 2,100 until 1,400 generations ago, then exponential growth to 500,000 |
| `himba` | Recent bottleneck: Ne falls from about 11,000 to about 500 between 86 and 12 generations ago, then grows to 1,000 |
| `expon` | Constant Ne until 150 generations ago, then exponential growth to 1,000,000 |
| `ooa2` | Two-population Out-of-Africa model (EUR and AFR, with migration) |
| `ashkenazi` | Ashkenazi Jewish model from stdpopsim (`AshkSub_7G19`) |

## Supplying your own tree sequences

Instead of `custom_demo`, set `custom_sim` to a function that returns a tree sequence for one chromosome:

```yaml
custom_sim:
  path: my_loader.py
  object: load_chromosome
```

```python
def load_chromosome(chrom, args):
    """chrom: chromosome number (1..end_chr); args: the run's YAML settings as a dict."""
    ...
    return ts, rate   # tskit.TreeSequence, recombination rate (float or msprime.RateMap)
```

Mutations are added by the pipeline. `ibd_sims/load_quebec.py` is a working example. No "true Ne" can be plotted for custom simulations.

## Writing a custom analysis

Subclass `PostProcessor` from `ibd_sims/post_process.py`:

```python
from post_process import PostProcessor

class MyAnalysis(PostProcessor):
    sub_config_key = "my_analysis"   # must match the YAML block name
    resource_fields = ["local", "workers", "mem_gb", "time_min"]

    def execute(self, wait=True):
        self._execute_helper()       # creates the numbered output folder (self.out_dir)
        if self.single_iter:
            self._single_iter(self.iter_n)
        else:
            self._execute_loop(wait=wait)

    def _single_iter(self, iter_n):
        cfg = self._get_sub_config()
        prefix = f"{self.path}/iter{iter_n}"
        # Read from: {prefix}.ibd.gz, {prefix}.map, {prefix}.tmrca.gz
        # Write to:  self.out_dir
```

Then add it to your config:

```yaml
post_process: ibdne,my_analysis

my_analysis:
  object: MyAnalysis
  path: my_analysis.py
  my_param: 42
  workers: 4
  time_min: 60
```

Things to know:

- `self._get_sub_config()` returns your YAML block, with each value as an attribute (`cfg.my_param`).
- `self._get_resource(name)` looks in your block first, then falls back to the top-level default.
- `resource_fields` lists settings that don't affect results. They are ignored when deciding whether an existing numbered output folder can be reused.
- `self._execute_loop()` runs every replicate, either locally or on Slurm depending on `local`.
- Override `is_iter_complete(iter_n)` to let re-runs skip finished replicates.
- `self.single_iter` is true when only one replicate is being processed (this happens at the end of each replicate during `run.py simulate`).

**Analyses that combine all replicates** (like `PostProcessIBDSummary`) should put that step behind `if not self.single_iter`. In single-replicate mode, `execute()` returns right after `_single_iter`, so the combined output is only written when the full loop runs. To refresh it after a single-replicate run, run `run.py postprocess` on the run folder; replicates that are already done are skipped via `is_iter_complete`.

## Using your own SNP density and allele frequencies

With hap-ibd (the default), simulated VCFs are thinned to a realistic SNP density and minor allele frequency spectrum, described by `ukb_snps.pkl` (derived from UK Biobank genotypes). To build your own from PLINK files:

```bash
python -m ibd_sims.maf_buckets \
    --afreq-chr1 /path/to/chr1.afreq \
    --bim-chr1 /path/to/chr1.bim \
    --output my_snps.pkl
```

Then set `maf_pickle: my_snps.pkl` in `setup.yaml`.

## How a simulation runs

For each replicate and chromosome:

1. The demographic model is loaded from `custom_demo`.
2. If `pedigree_mode` is true, a Wright-Fisher pedigree covering the last `gen_end` generations is generated with the chosen mating model (`wf_pedigree.py`) and saved as `iter{n}_WF.pedigree`. msprime simulates through the pedigree, then switches to the coalescent for older generations. Otherwise the whole history is simulated with the coalescent.
3. Mutations are added at rate 1e-8.
4. IBD is called, either by writing a thinned VCF and running hap-ibd, or directly from the tree sequence (`tskit_ibd: true`, `ibd_from_ts.py`).
5. Each IBD segment is annotated with its TMRCA from the tree sequence.
6. Per-chromosome files are combined into `iter{n}.ibd.gz`, `iter{n}.map` and `iter{n}.tmrca.gz`, and intermediate files are deleted unless `keep_all_files: true`.

Random seeds are reproducible: each replicate's seed is derived from the run folder path and replicate number, and each chromosome gets its own seed from that.

## Repository layout

```
├── run.py                        # main entry point: simulate / postprocess
├── setup.yaml                    # paths to external tools on your machine
├── environment.yaml              # conda environment
├── ukb_snps.pkl                  # SNP density / allele frequency profile
├── yaml_files/                   # example simulation and experiment configs
├── docs/                         # detailed documentation
├── paper/                        # scripts for the paper's analyses
└── ibd_sims/                     # pipeline source code
    ├── simulate.py               # job orchestration (local or Slurm)
    ├── simulations.py            # core simulation: msprime, VCF, hap-ibd, TMRCA
    ├── ibd_from_ts.py            # IBD straight from the tree sequence (tskit_ibd)
    ├── demography.py             # built-in demographic models
    ├── wf_pedigree.py            # Wright-Fisher pedigree generation
    ├── write_vcf.py              # VCF and genetic map output with SNP thinning
    ├── maf_buckets.py            # build a SNP density / allele frequency profile
    ├── post_process.py           # analysis orchestration, PostProcessor base class
    ├── post_modules.py           # built-in analyses: IBDNe, HapNe, IBD summary, purple nodes
    ├── run_hapne.py              # HapNe-IBD and HapNe-LD helpers
    ├── filter_ibd.py             # related / unrelated / random sample subsets
    ├── purple.py                 # purple-node matrix
    ├── concat_tmrca.py           # combine per-chromosome TMRCA files
    ├── experiment.py             # simulation experiment manager
    ├── postprocess_experiment.py # post-processing experiment manager
    ├── analyze_experiment.py     # load Ne estimates across an experiment
    ├── plot_Ne.py                # plot Ne estimates against the truth
    └── utils.py                  # --set override handling
```
