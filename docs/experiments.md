# Running many simulations: experiments

When you want to compare several demographic histories, mating models or genome layouts, writing one YAML per combination gets tedious. The two experiment managers do it for you:

1. **Simulation experiments** (`ibd_sims/experiment.py`): describe the combinations once, and it writes one simulation config per combination and tracks their progress.
2. **Post-processing experiments** (`ibd_sims/postprocess_experiment.py`): run the same set of analyses (with every combination of settings you list) across all simulations in an experiment.

Both work the same way: **describe** (preview, writes nothing) → **init** (create files) → **commands** (print the commands to run) → **status** (check progress).

All commands are run from the repository root.

## Simulation experiments

### 1. Write an experiment file

See `yaml_files/experiment.yaml` for a working example.

```yaml
experiment: my_experiment      # output folder name

# Shared by every simulation
iter: 50
samples: 1000
sim_workers: 30

# Default resources (can be overridden below)
gb: 8
sim_min: 30
nthreads: 8
keep_all_files: false

# ── Things to vary: every combination is simulated ──
end_chr:
  22: {}                       # no resource overrides
  30:
    resources:
      sim_min: 45              # override for this genome layout

demographies:
  constant_Ne_10k:
    object: constant_Ne
    path: ibd_sims/demography.py

  constant_Ne_100k:
    object: constant_Ne100k
    path: ibd_sims/demography.py
    resources:
      sim_min: 45

mating:
  DTWF_di:
    pedigree_mode: true
    mating: di
    gen_end: 25
    pedigree_file: null

# Custom simulations are added once each (not combined with end_chr or mating)
custom_sims:
  quebec:
    object: load_random_10000
    path: ibd_sims/load_quebec.py
    resources:
      gb: 32
      end_chr: 22

post_processing: ibdne,hapne_ibd
```

When several of these set the same resource (say `sim_min`), the largest value wins.

### 2. Preview the plan

```bash
python ibd_sims/experiment.py describe yaml_files/experiment.yaml
```

### 3. Create the configs

```bash
python ibd_sims/experiment.py init yaml_files/experiment.yaml
```

This writes one simulation config per combination:

```
my_experiment/
└── yaml_files/
    ├── constant_Ne_10k__DTWF_di__chr22.yaml
    ├── constant_Ne_10k__DTWF_di__chr30.yaml
    ├── constant_Ne_100k__DTWF_di__chr22.yaml
    ├── constant_Ne_100k__DTWF_di__chr30.yaml
    └── quebec.yaml
```

### 4. Get the commands to run

```bash
# All simulations, submitted to Slurm in parallel (--no-wait is the default here)
python ibd_sims/experiment.py commands yaml_files/experiment.yaml

# Only simulations that haven't finished
python ibd_sims/experiment.py commands yaml_files/experiment.yaml --pending-only

# One at a time: wait for each simulation to finish before starting the next
python ibd_sims/experiment.py commands yaml_files/experiment.yaml --wait
```

This prints one `python run.py simulate ...` command per simulation. Copy them into your terminal or a script.

### 5. Check progress

```bash
python ibd_sims/experiment.py status yaml_files/experiment.yaml
```

### Adding analyses afterwards

Either edit the generated configs in `my_experiment/yaml_files/` and run `python run.py postprocess my_experiment/<simulation>/`, or use a post-processing experiment (below) to do all simulations at once.

## Post-processing experiments

### 1. Write a post-processing file

See `yaml_files/postprocess_experiment.yaml` for a working example.

```yaml
experiment_directory: my_experiment

postprocess: [ibdne, hapne_ibd, ibd_summary]

ibdne:
  path: ibd_sims/post_modules.py
  object: PostProcessIBDNe
  mincm: 2
  trimcm: 0.2
  gmin: 1
  gmax: 300
  nboots: 80
  nits: 1000
  npairs: 0
  workers: 8
  mem_gb: 16
  time_min: 120

  # Every combination of these is run
  combo_args:
    filtersamples: [true, false]
    filter: [null, related, unrelated]

  # Add or remove specific combinations by hand
  add_combo:
    combo1:
      filtersamples: true
      filter: null

  ignore_combo:
    combo1:
      filtersamples: true
      filter: null

hapne_ibd:
  path: ibd_sims/post_modules.py
  object: PostProcessHapNeIBD
  workers: 4
  mem_gb: 16
  time_min: 120
  combo_args:
    filter: [null, related, unrelated]

ibd_summary:
  path: ibd_sims/post_modules.py
  object: PostProcessIBDSummary
  mem_gb: 1
  time_min: 10
  workers: 1
  filters: [null, related, unrelated]
```

### 2. Preview the plan

```bash
python ibd_sims/postprocess_experiment.py describe yaml_files/postprocess_experiment.yaml
```

### 3. Initialise

```bash
python ibd_sims/postprocess_experiment.py init yaml_files/postprocess_experiment.yaml
```

This creates two files in `my_experiment/`:

- `postprocess.tsv`: one row per (analysis, settings combination), with a `status` column (`new`, `rerun` or `complete`)
- `postprocess.yaml`: the shared settings (everything except the combination axes), used by the generated commands

### 4. Get the commands to run

```bash
# Run locally (one command after another)
python ibd_sims/postprocess_experiment.py commands yaml_files/postprocess_experiment.yaml

# Submit to Slurm instead
python ibd_sims/postprocess_experiment.py commands yaml_files/postprocess_experiment.yaml --no-local

# Submit to Slurm without waiting for each command to finish
python ibd_sims/postprocess_experiment.py commands yaml_files/postprocess_experiment.yaml --no-local --no-wait
```

This prints one `python run.py postprocess ...` command per (simulation, analysis combination) and also writes them to `my_experiment/postprocess_scripts/run.sh`. Only rows with status `new` or `rerun` are included.

### 5. Check progress

```bash
python ibd_sims/postprocess_experiment.py status yaml_files/postprocess_experiment.yaml
```

This checks the output files and prints a table:

```
postprocess  directory      progress  status
---------------------------------------------
ibdne        ibdne/001      45/50     rerun
ibdne        ibdne/002      50/50     complete
hapne_ibd    hapne_ibd/001  0/50      new
```

It also updates `postprocess.tsv`, so running `commands` again picks up anything unfinished.

## Plotting an experiment

```bash
python ibd_sims/plot_Ne.py my_experiment/
```

This writes one plot per simulation to `my_experiment/plots/<simulation>_Ne_plot.png`, showing every IBDNe and HapNe-IBD estimate against the true Ne. Lines are coloured by tool (greens for IBDNe, oranges for HapNe-IBD) and labelled with the analysis settings from `postprocess.tsv`. With more than 10 replicates, a 5th–95th percentile band is drawn. Replicates whose largest Ne estimate is over 100× the median are left out as outliers.

| Option | What it does |
|--------|--------------|
| `--no-vlines` | Hide the log2(Ne) reference lines. |
| `--save-pickle` | Also save the plotted data to `plots/Ne_data.pkl`, for re-plotting in a notebook. |
