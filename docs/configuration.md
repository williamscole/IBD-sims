# Configuration reference

Every run is driven by one YAML file. This page lists every setting. For a gentler introduction, see the [README](../README.md).

- [Command-line options](#command-line-options)
- [Simulation settings](#simulation-settings)
  - [Choosing a genome (`end_chr`)](#choosing-a-genome-end_chr)
  - [`tskit_ibd`: calling IBD from the tree sequence](#tskit_ibd-calling-ibd-from-the-tree-sequence)
  - [`sim_workers` and job modes](#sim_workers-and-job-modes)
  - [Slurm queue limits (`--max-jobs`)](#slurm-queue-limits---max-jobs)
- [Post-processing settings](#post-processing-settings)
  - [Sample filtering](#sample-filtering)
  - [Built-in analyses](#built-in-analyses)
- [Machine settings (`setup.yaml`)](#machine-settings-setupyaml)

## Command-line options

All commands are run from the repository root.

### `python run.py simulate CONFIG_OR_RUN_DIR`

| Option | What it does |
|--------|--------------|
| `--local` | Run on this machine instead of submitting to Slurm. |
| `--workers N` | With `--local`, run at most `N` jobs at once (default: all CPUs). |
| `--set KEY=VALUE ...` | Override settings from the YAML, e.g. `--set iter=5 pedigree.mating=mono`. Use dots for nested keys. Can be repeated. Lists can't be set this way. |
| `--no-wait` | Submit the Slurm jobs and exit straight away. Post-processing is skipped; run `run.py postprocess` later. |
| `--max-jobs N` | Maximum number of Slurm jobs to queue at once (default 1000). See [below](#slurm-queue-limits---max-jobs). |

Passing an existing run folder instead of a YAML file resumes that run, skipping chromosomes that already finished.

### `python run.py postprocess RUN_DIR`

| Option | What it does |
|--------|--------------|
| `--set KEY=VALUE ...` | Override analysis settings, e.g. `--set ibdne.nboots=100 ibdne.mincm=3`. |
| `--set local=false` | Submit the analyses to Slurm. **Post-processing runs locally by default.** |
| `--no-wait` | With Slurm, submit jobs and exit straight away. |

## Simulation settings

```yaml
# ── Output ───────────────────────────────────────────────
base_dir: my_runs           # parent folder for output (null = current folder)
label: constant_Ne_n1000    # run folder name; a multi-line label is joined with "_".
                            # If the folder exists, _001, _002, ... is appended.
keep_all_files: false       # also keep the VCF and HBD files (needs bcftools)
keep_trees: false           # also save each chromosome's tree sequence (.trees)
tskit_ibd: false            # true = call IBD from the tree sequence (see below)

# ── Computational resources (simulation) ─────────────────
gb: 8                       # memory (GB) per simulation job
sim_min: 30                 # time limit (minutes) per chromosome
nthreads: 8                 # threads for hap-ibd (unused with tskit_ibd)
sim_workers: 1              # chromosomes simulated in parallel inside one job

# ── What to simulate ─────────────────────────────────────
iter: 50                    # number of independent replicates
samples: 1000               # number of diploid individuals
end_chr: 30                 # genome layout; see "Choosing a genome" below
subsample_frac: 0.25        # size of the random/related/unrelated subsets,
                            # as a fraction of samples (0.25 = N/4)

# ── Demographic history (set one of these two) ──────────
custom_demo:
  path: ibd_sims/demography.py   # Python file with an msprime.Demography object
  object: constant_Ne            # name of that object
custom_sim:
  path: null                     # Python file with a function returning (ts, rate)
  object: null                   # name of that function

# ── Mating model ─────────────────────────────────────────
pedigree:
  pedigree_mode: true       # true = explicit Wright-Fisher pedigree for recent generations
  mating: di                # "di" = random mating, "mono" = monogamous
  gen_end: 25               # pedigree depth (generations) before switching to the coalescent
  pedigree_file: null       # use an existing pedigree file instead of generating one
```

Older configs may contain `dir_name`; `run.py simulate` ignores it (the folder name comes from `label`).

### Choosing a genome (`end_chr`)

| `end_chr` | Genome | Recombination |
|-----------|--------|---------------|
| `1` | one 50 Mb chromosome | constant 1e-8 per bp |
| `2` | two 50 Mb chromosomes | constant 1e-8 per bp |
| `30` | thirty 100 Mb chromosomes | constant 1e-8 per bp |
| `22` | the 22 human autosomes | HapMap GRCh37 maps (needs `hapmap_chr1` in `setup.yaml`) |

Other values are not supported.

### `tskit_ibd`: calling IBD from the tree sequence

By default, each chromosome is written to a thinned VCF and IBD is detected with hap-ibd. Setting `tskit_ibd: true` skips both steps and calls IBD directly from the simulated tree sequence, avoiding the cost of writing a VCF and running hap-ibd. The output has the same format as hap-ibd's `.ibd.gz`, so post-processing needs no changes.

How it differs from the hap-ibd route:

- **Exact IBD.** Segments come from the true genealogy, so there is no SNP thinning, phasing, or detection error. Use the default route when you want IBD with realistic detection error.
- **Segment boundaries** are the exact tree-sequence breakpoints rather than the first and last SNP of a segment.
- **Segments** are stretches with a single most recent common ancestor, kept if they are at least 2 cM long (`min_cm`); touching segments are merged.
- **No VCF or HBD output.** `keep_all_files` has no effect on these, and bcftools is not required.
- **Genetic map.** With no SNPs, the `.map` file is built from the recombination rate: a constant rate gives rows every 10 kb along the chromosome, and an `msprime.RateMap` gives the map's own breakpoints (merged with the same 10 kb grid). Every called segment's start and end is also added as a row, and the last row is always the end of the chromosome.
- **Not used:** `hap_ibd_jar` and `maf_pickle` in `setup.yaml`, and `nthreads`.
- **TMRCA annotations** are computed as in the default route.
- **Genotype-based analyses are unavailable** (HapNe-LD needs genotypes). IBD-based post-processing (IBDNe, HapNe-IBD, `purple_nodes`, `ibd_summary`) is unaffected.

### `sim_workers` and job modes

The pipeline picks one of two ways to split up the work:

- **Per-iteration mode** (`iter >= 3`): one job per replicate, running all chromosomes. `sim_workers` sets how many chromosomes run in parallel inside that job; on Slurm the job asks for `sim_workers` CPUs and `gb × sim_workers` memory.
- **Per-chromosome mode** (`iter < 3`): one job per (replicate, chromosome) pair. `sim_workers` is not used.

### Slurm queue limits (`--max-jobs`)

In per-chromosome mode, if the number of jobs exceeds `--max-jobs` (default 1000), jobs are submitted in batches of about `max_jobs / 4`, and each batch waits for the previous one to finish. `--no-wait` is ignored when batching is needed.

## Post-processing settings

`post_process` is a comma-separated list of analyses to run, or `null` for none. Each analysis has its own block, named after it, with:

- `path` and `object`: the Python file and class that implement it (use the values shown below for built-in analyses)
- its analysis settings
- optional resources: `workers`, `mem_gb`, `time_min`. Anything not set falls back to top-level `workers`, `mem_gb`, `time_min` and `local`.

Run analyses with `run.py postprocess RUN_DIR`, any time after simulating. Analyses listed in `post_process` also run automatically at the end of `run.py simulate`, but currently only when `iter` is 1 or 2 (per-chromosome mode).

Each run of an analysis writes to a numbered folder (`ibdne/001/`, `ibdne/002/`, ...). Re-running with the same settings reuses the matching folder, skipping replicates that already finished. Changing any non-resource setting creates a new folder. Each folder's `args.yaml` records the settings used.

### Sample filtering

The `filter` setting (IBDNe, HapNe-IBD, HapNe-LD, IBD summary) controls which individuals are included:

| Value | Who is included |
|-------|-----------------|
| `null` / `none` | Everyone |
| `random` | A random subset of `round(samples × subsample_frac)` people |
| `related` | A subset of the same size, enriched for 1st–3rd degree relative pairs |
| `unrelated` | A subset of the same size, with 1st–3rd degree relatives pruned out |

The chosen individuals are cached next to the IBD file as `iter{n}_random.txt`, `iter{n}_related.txt` and `iter{n}_unrelated.txt`, and reused on later runs.

### Built-in analyses

**IBDNe** (`ibdne`): estimates Ne over time. Needs Java and `ibdne_jar`.

```yaml
ibdne:
  path: post_modules.py
  object: PostProcessIBDNe
  filter: null              # sample filtering (see above)
  filtersamples: false      # IBDNe's own filtersamples option
  mincm: 2                  # minimum segment length (cM)
  trimcm: 0.2               # trim this much off each segment end (cM)
  gmin: 1                   # first generation to estimate
  gmax: 300                 # last generation to estimate
  nboots: 80                # bootstrap replicates
  nits: 1000                # iterations
  npairs: 0                 # max pairs (0 = all)
  workers: 8
  mem_gb: 16
  time_min: 120
```

**HapNe-IBD** (`hapne_ibd`): estimates Ne over time with HapNe. Needs the [HapNe patches](../README.md#optional-hapne).

```yaml
hapne_ibd:
  path: post_modules.py
  object: PostProcessHapNeIBD
  filter: null
  workers: 4
  mem_gb: 16
  time_min: 120
```

**HapNe-LD** (`hapne_ld`): estimates Ne from linkage disequilibrium. Needs genotypes (so not `tskit_ibd: true`) and the HapNe patches. Currently slow and may not work.

```yaml
hapne_ld:
  path: post_modules.py
  object: PostProcessHapNeLD
  filter: null
  workers: 4
  mem_gb: 16
  time_min: 120
```

**IBD summary** (`ibd_summary`): for each replicate and filter, counts samples, segments and sharing pairs, and reports mean pairwise and total IBD (cM). Needs no external tools. Each replicate writes `iter{n}.tsv`; `ibd_summary.tsv` combines all replicates. The combined file is only written by `run.py postprocess`, not by the per-replicate step at the end of `run.py simulate`.

```yaml
ibd_summary:
  path: post_modules.py
  object: PostProcessIBDSummary
  mincm: 2                              # minimum segment length (cM)
  filters: [null, related, unrelated]   # filters to summarise
  workers: 1
  mem_gb: 1
  time_min: 10
```

**Purple nodes** (`purple_nodes`): computes the purple-node matrix for each replicate (`iter{n}.npy`).

```yaml
purple_nodes:
  path: post_modules.py
  object: PostProcessPurple
  workers: 4
  mem_gb: 8
  time_min: 60
```

## Machine settings (`setup.yaml`)

`setup.yaml` in the repository root holds paths specific to your machine. Only fill in what you use.

| Key | Needed when |
|-----|-------------|
| `maf_pickle` | Using hap-ibd (the default). `ukb_snps.pkl` is included. |
| `hap_ibd_jar` | Using hap-ibd (the default). |
| `ibdne_jar` | Running IBDNe. |
| `hapmap_chr1` | `end_chr: 22`. Path to the chromosome 1 HapMap-format genetic map (GRCh37); other chromosomes are found by replacing `chr1` with `chr{n}`. |
| `quebec_ts_chr1` | Running the Quebec config (`arg7.yaml`). |
