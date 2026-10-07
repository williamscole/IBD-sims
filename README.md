# IBD-sims

Simulate realistic **identity-by-descent (IBD) segments** under a population history you choose, then test how well tools like [IBDNe](https://faculty.washington.edu/browning/ibdne.html) recover that history.

You pick a demographic history (constant size, Out-of-Africa, a bottleneck, or your own) and a mating model (random or monogamous). The pipeline then:

1. **Simulates** genomes for a sample of people and finds the IBD segments they share.
2. **Analyses** those segments, for example estimating effective population size (Ne) over time with IBDNe or HapNe.
3. **Plots** the estimates against the true history you simulated.

Because the truth is known, you can see exactly where an Ne-estimation method gets things right and wrong. The IBD segments are plain text files, so you can also use them for any other analysis.

<details>
<summary>What's under the hood?</summary>

Genealogies are simulated with [msprime](https://tskit.dev/msprime/), either with the standard coalescent or with an explicit Wright-Fisher pedigree for the most recent generations. IBD segments are found either with [hap-ibd](https://github.com/browning-lab/hap-ibd) on simulated genotypes (realistic detection error), or read exactly from the simulated genealogy with tskit (fast, no external tools). Jobs run on your own machine or on a Slurm cluster.
</details>

## Contents

- [Quick start (about 5 minutes)](#quick-start-about-5-minutes)
- [Full installation](#full-installation)
- [Running simulations](#running-simulations)
- [Estimating Ne and other analyses](#estimating-ne-and-other-analyses)
- [Plotting](#plotting)
- [Understanding the output](#understanding-the-output)
- [Writing your own config](#writing-your-own-config)
- [Troubleshooting](#troubleshooting)
- [More documentation](#more-documentation)
- [Known limitations and to-do](#known-limitations-and-to-do)

## Quick start (about 5 minutes)

This runs a tiny simulation on your own computer. It needs only [conda](https://docs.conda.io/en/latest/miniconda.html): no Java, no cluster, no extra downloads.

**1. Get the code and install the Python environment** (the environment step takes a while the first time):

```bash
git clone https://github.com/williamscole/IBD-sims.git
cd IBD-sims
conda env create -f environment.yaml
conda activate ibd-sims
```

**2. Let Python find the pipeline code** (one-time step):

```bash
echo "$(cd ibd_sims && pwd)" > $(python -c "import site; print(site.getsitepackages()[0])")/ibd-sims.pth
```

**3. Run the example:**

```bash
python run.py simulate yaml_files/quickstart.yaml --local --workers 2
```

This simulates 200 people on two chromosomes and summarises the IBD they share. When it prints `Done.`, look in `quickstart/quickstart/`:

- `iter1.ibd.gz`: the IBD segments (one per line)
- `ibd_summary/001/iter1.tsv`: a small table of how much IBD was found

That's it, the pipeline works. Next, follow the [full installation](#full-installation) to set up the external tools for IBD calling and Ne estimation.

## Full installation

### 1. Python environment

Same as the quick start: create the conda environment and register the pipeline (steps 1 and 2 above).

> **Tip:** run every command from the repository root (the `IBD-sims/` folder). Config files refer to paths like `ibd_sims/demography.py` relative to it.

### 2. External tools: install only what you need

| Tool | Needed for | Not needed if |
|------|------------|---------------|
| Java 8+ | hap-ibd and IBDNe | you use neither |
| [hap-ibd](https://github.com/browning-lab/hap-ibd) jar | Calling IBD from simulated genotypes (the default) | you set `tskit_ibd: true` |
| [IBDNe](https://faculty.washington.edu/browning/ibdne.html) jar | Estimating Ne with IBDNe | you don't run IBDNe |
| bcftools | Keeping the combined VCF (`keep_all_files: true`) | you leave `keep_all_files: false` |
| HapMap GRCh37 genetic maps | Simulating the real human autosomes (`end_chr: 22`) | you use simulated chromosomes (`end_chr: 30`, the default in the examples) |

### 3. Tell the pipeline where the tools are

Edit `setup.yaml` in the repository root. Leave entries you don't need as they are.

```yaml
maf_pickle: ukb_snps.pkl                                # included; leave as is
hap_ibd_jar: /path/to/hap-ibd.jar
ibdne_jar: /path/to/ibdne.jar
hapmap_chr1: /path/to/genetic_map_GRCh37_chr1.txt.gz    # only for end_chr: 22
```

### Optional: HapNe

HapNe is installed by the conda environment, but the pinned version needs two small fixes to work with recent NumPy. You only need these to run HapNe-IBD or HapNe-LD.

```bash
HAPNE=$(python -c "import hapne, os; print(os.path.dirname(hapne.__file__))")
sed -i 's/times\[ii + 1\] = t_quantile/times[ii + 1] = np.asarray(t_quantile).item()/' $HAPNE/backend/DemographicHistory.py
sed -i 's/n = n.ravel()/n = np.asarray(n).ravel()/' $HAPNE/utils.py
```

## Running simulations

Everything goes through `run.py`. A simulation is described by a YAML config file; `yaml_files/` has ready-made ones.

```bash
# Run on a Slurm cluster (the default) and wait for it to finish
python run.py simulate yaml_files/arg1.yaml

# Run on this computer, at most 8 jobs at a time
python run.py simulate yaml_files/arg1.yaml --local --workers 8

# Submit to Slurm and return immediately
python run.py simulate yaml_files/arg1.yaml --no-wait
```

**Change settings without editing the file** with `--set`:

```bash
python run.py simulate yaml_files/arg1.yaml --set iter=5 pedigree.mating=mono
```

**Pick up where you left off.** If a run was interrupted, pass its output folder instead of a config. Finished chromosomes are skipped.

```bash
python run.py simulate path/to/run_folder/
```

Each run gets its own output folder, named after the config's `label`. Running the same config again creates a new folder (`..._001`, `..._002`) rather than overwriting.

See [all command-line options](docs/configuration.md#command-line-options), including `--max-jobs` for clusters with queue limits.

### Ready-made configs

| Config | History | People | Mating model | Replicates |
|--------|---------|--------|--------------|------------|
| `quickstart.yaml` | Constant Ne 10,000 | 200 | Coalescent | 1 (tiny, no external tools) |
| `debug.yaml` | Constant Ne 10,000 | 1,000 | Coalescent | 1 |
| `arg1.yaml` | Constant Ne 10,000 | 1,000 | Random, 25-generation pedigree | 50 |
| `arg2.yaml` | Constant Ne 10,000 | 1,000 | Monogamous, 25-generation pedigree | 50 |
| `arg3.yaml` | Constant Ne 100,000 | 1,000 | Random, 25-generation pedigree | 50 |
| `arg4.yaml` | Constant Ne 100,000 | 1,000 | Monogamous, 25-generation pedigree | 50 |
| `arg5.yaml` | Out-of-Africa (2 populations) | 2,000 | Random, 25-generation pedigree | 50 |
| `arg6.yaml` | Out-of-Africa (2 populations) | 2,000 | Monogamous, 25-generation pedigree | 50 |
| `arg7.yaml` | Quebec (empirical tree sequences) | 10,000 | n/a | 1 |
| `arg8.yaml` | Ashkenazi | 1,000 | Random, 12-generation pedigree | 50 |

Apart from `quickstart.yaml`, these call IBD with hap-ibd and run no analyses by default. Add `--set tskit_ibd=true` to skip hap-ibd. The other files in `yaml_files/` are older configs or [experiment](docs/experiments.md) files.

## Estimating Ne and other analyses

After simulating, you run analyses ("post-processing") on the IBD segments. Choose them with `post_process` in the config:

```bash
python run.py postprocess path/to/run_folder/ --set post_process=ibdne
```

You can re-run analyses as often as you like with different settings, without re-simulating.

| Analysis | What it does | Needs |
|----------|--------------|-------|
| `ibdne` | Estimates Ne over time with IBDNe | Java, `ibdne_jar` |
| `hapne_ibd` | Estimates Ne over time with HapNe-IBD | [HapNe fixes](#optional-hapne) |
| `hapne_ld` | Estimates Ne from linkage disequilibrium (experimental) | HapNe fixes, genotypes (not `tskit_ibd`) |
| `ibd_summary` | Counts segments and sharing pairs; total and mean IBD | nothing |
| `purple_nodes` | Computes the purple-node matrix | nothing |

Run several at once with a comma-separated list: `--set post_process=ibdne,ibd_summary`.

**Tweaking settings:** each analysis has a block in the config (e.g. `ibdne:`) with its settings. Override them like this:

```bash
python run.py postprocess path/to/run_folder/ --set post_process=ibdne ibdne.mincm=3 ibdne.nboots=100
```

**Numbered results:** results go in numbered folders, e.g. `ibdne/001/`, `ibdne/002/`. Running again with the same settings reuses the same folder (finished replicates are skipped); changing a setting makes a new one. Each folder has an `args.yaml` recording the exact settings used.

**Subsets of people:** most analyses take a `filter` setting to analyse everyone (`null`), or a random, related-enriched or unrelated subset.

**Where it runs:** post-processing runs on your computer by default, even if the simulation ran on Slurm. Add `--set local=false` to submit it to Slurm.

Every setting for every analysis is in the [configuration reference](docs/configuration.md#post-processing-settings).

## Plotting

Plotting works on [experiments](docs/experiments.md) (a folder of related runs):

```bash
python ibd_sims/plot_Ne.py my_experiment/
```

This saves one plot per simulation in `my_experiment/plots/`, comparing every IBDNe and HapNe-IBD estimate with the true Ne. See [plotting an experiment](docs/experiments.md#plotting-an-experiment) for options.

## Understanding the output

A finished run folder looks like this (one set of `iter` files per replicate):

```
quickstart/quickstart/
├── args.yaml            # the exact settings used for this run
├── iter1.ibd.gz         # IBD segments
├── iter1.map            # genetic map (PLINK format: chr, id, cM, bp)
├── iter1.tmrca.gz       # when each segment's common ancestor lived (TMRCA)
├── iter1_related.txt    # people chosen for the "related" subset (also _random, _unrelated)
├── ibdne/001/           # results of each analysis run, in numbered folders
├── ibd_summary/001/
├── slurm/               # job logs: look here first if something fails
└── errors/
```

Each line of `iter{n}.ibd.gz` is one segment shared by two people, in [hap-ibd's format](https://github.com/browning-lab/hap-ibd#output-files):

```
id1  hap1  id2  hap2  chromosome  start_bp  end_bp  length_cM
```

## Writing your own config

The easiest start is to copy `yaml_files/quickstart.yaml` (small and commented) or one of the `arg*.yaml` files. These are the settings you'll change most:

| Setting | What it controls | Example |
|---------|------------------|---------|
| `label` / `base_dir` | Output folder name and parent folder | `label: my_run` |
| `iter` | Number of independent replicates | `50` |
| `samples` | Number of people (diploid individuals) | `1000` |
| `end_chr` | Genome: `30` = thirty 100 Mb chromosomes, `22` = human autosomes, `1`/`2` = one or two 50 Mb chromosomes | `30` |
| `custom_demo` | Demographic history: a file and an `msprime.Demography` in it | `object: ooa2` |
| `pedigree.pedigree_mode` | Use an explicit pedigree for recent generations | `true` |
| `pedigree.mating` | `di` (random) or `mono` (monogamous) | `di` |
| `tskit_ibd` | `true` = exact IBD from the genealogy (fast, no hap-ibd); `false` = hap-ibd on genotypes (realistic errors) | `false` |
| `gb`, `sim_min` | Memory (GB) and time limit (minutes) per chromosome | `8`, `30` |
| `post_process` | Analyses to run after simulating | `ibdne,ibd_summary` |

Built-in histories (in `ibd_sims/demography.py`): `constant_Ne`, `constant_Ne100k`, `euro_bottleneck`, `himba`, `expon`, `ooa2`, `ashkenazi`. You can also [add your own](docs/extending.md#adding-a-demographic-model).

Every setting is in the [configuration reference](docs/configuration.md).

## Troubleshooting

**`conda env create` fails (e.g. on macOS):** `environment.yaml` pins exact Linux builds. Create a plain environment instead:

```bash
conda create -n ibd-sims python=3.12
conda activate ibd-sims
pip install msprime tskit stdpopsim submitit pyyaml numpy pandas scipy matplotlib seaborn networkx tszip polars
pip install hapne==1.20240807 pandas-plink   # only needed for HapNe
```

**A job failed. Where do I look?** Open `<run_folder>/slurm/*_log.err` (each job's error log, used for local runs too) and `<run_folder>/errors/`.

**`ModuleNotFoundError: No module named 'simulations'`** (or `post_modules`, `simulate`, ...): Python can't find the pipeline code. Run step 2 of the [quick start](#quick-start-about-5-minutes) inside the activated `ibd-sims` environment.

**`FileNotFoundError` for `ibd_sims/demography.py` or `ukb_snps.pkl`:** run commands from the repository root.

**`UnboundLocalError: ... 'sequence_length'`:** `end_chr` must be 1, 2, 22 or 30.

**hap-ibd or IBDNe errors:** check that the jar paths in `setup.yaml` are right and that `java -version` works. To skip hap-ibd entirely, use `--set tskit_ibd=true`.

**`bcftools not found`:** only needed with `keep_all_files: true`. Install bcftools or set it to `false`.

**Jobs run out of memory:** increase `gb` (simulation) or the analysis's `mem_gb`.

**Too many Slurm jobs for my cluster's queue limit:** use `--max-jobs`, e.g. `--max-jobs 200`.

## More documentation

- [Configuration reference](docs/configuration.md): every setting and command-line option
- [Experiments](docs/experiments.md): running and plotting many simulations at once
- [Extending the pipeline](docs/extending.md): your own demographic models, tree sequences and analyses; how a simulation runs; repository layout
- [`llm.txt`](llm.txt): a single-file summary to paste into an AI assistant when you need help

## Known limitations and to-do

- `python run.py plot` is out of date and fails; use `python ibd_sims/plot_Ne.py` on an experiment folder instead.
- HapNe-LD is slow and may not work.
- `tskit_ibd: true` is not chosen automatically when no hap-ibd jar is configured.
- `ibd_sims/experiment.py commands` has no global `--max-jobs`, so with `--no-wait` (its default) a large experiment can exceed a cluster's per-user Slurm queue limit.
- Resuming re-runs every unfinished replicate; you can't yet pick specific replicates.
- Analyses listed in a config only run automatically during `run.py simulate` when `iter` is 1 or 2. With more replicates, run `run.py postprocess` afterwards.
- Custom simulations (`custom_sim`) set `end_chr` through a workaround (under `resources` in experiment files).
- Long-term: integrate ped-sim for more realistic IBD between close relatives.
