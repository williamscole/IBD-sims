# min_span sensitivity experiment: how it was run

This experiment asks how much the IBD segments called directly from the tree sequence depend on
`min_span`, and how they compare with hap-ibd calls made on the **same** simulated tree sequences.
It compares segment-level agreement (RMSE and bias of per-pair total IBD, segment counts) and the
downstream IBDNe Ne estimate.

Everything is run from the repository root. Paths below are relative to it.

## Design

| | |
|---|---|
| Replicates | 25 (`iter: 25`) |
| Sample | 100 diploids (200 haplotypes) |
| Genome | 30 chromosomes of 100 Mb, constant recombination rate 1e-8 (`end_chr: 30`) |
| Demography | `constant_Ne` from `ibd_sims/demography.py` (Ne = 10,000), coalescent (`pedigree_mode: false`) |
| Reference | tskit calls at `min_span` = 20 kb |
| tskit settings compared | `min_span` = 20, 50, 100, 200, 500, 1000 kb |
| Other caller | hap-ibd on a VCF written from the same tree sequence |
| Minimum segment length | 2 cM (all callers) |

**Same tree for every caller.** The simulation seed is a deterministic function of the run
directory, iteration and chromosome (`sim()` in `ibd_sims/simulations.py`). Each tree sequence is
simulated once. hap-ibd's VCF is written from it, and every tskit `min_span` is then called from the
saved `.trees` file. So the callers differ only in how IBD is called, never in the data.

**tskit segment definition.** `ts.ibd_segments(min_span=...)` returns pieces of constant MRCA per
sample-node pair. `ibd_from_ts.py` merges adjacent pieces with the same MRCA. `min_span` is applied
per piece before merging, so a large `min_span` undercounts segments built from many small pieces.
That is the effect this experiment measures.

## Prerequisites

1. `setup.yaml` has the right `hap_ibd_jar`, `ibdne_jar`, `hapmap_chr1` and `maf_pickle`.
2. `ibd_sims/ibd_from_ts.py` is the fixed version: `ibd_from_ts(...)` takes `same_mrca` and
   `min_span`, and the default is not `min_span=0`. `minspan_sweep.py call` warns if `same_mrca` is
   missing.
3. `sim()` in `ibd_sims/simulations.py` saves the tree sequence when `keep_trees: true`, i.e. the
   cluster copy calls `add_tmrca(prefix, ts, yargs.get("keep_trees", False))` so that
   `iter{n}_chr{c}.trees` is written. (The copy of the repo this README was written against does not
   have this change, since it passes `False`.) Without it there are no `.trees` files and `call` has
   nothing to read.
4. `sim()` times the hap-ibd route so the time is in the log:

   ```python
   else:
       t1 = time.time()
       write_vcf(ts, prefix, chrom, rate, seed, snps_pkl=config["maf_pickle"])
       run_hapibd(prefix, yargs["gb"], hapibd_jar=config["hap_ibd_jar"])
       print(f"Time to write VCF and call IBD: {round(time.time() - t1, 4)} (iter={iter_n}, chrom={chrom})")
   ```

   The `print` must be inside the `else:` and `import time` must be at the top of the file.
5. Python packages from `environment.yaml` (msprime, tskit, numpy, pandas, matplotlib, PyYAML,
   submitit), `java` for hap-ibd and IBDNe, and Slurm.

## Step 0: start clean

The rerun must regenerate everything so the timing lines are in the logs. Move aside any earlier
output rather than mixing it with the new run:

```bash
mv minspan_sweep minspan_sweep_old      # or delete it
```

## Step 1: the experiment definition

`yaml_files/minspan_sweep.yaml`:

```yaml
experiment: minspan_sweep

# Fixed across all simulations
iter: 25
samples: 100
sim_workers: 10

# Resources (per chromosome; the per-iteration job multiplies these)
gb: 8
sim_min: 30
nthreads: 2

# This experiment
tskit_ibd: false        # hap-ibd route: gives the hap-ibd .ibd.gz to compare against
keep_trees: true        # save iter{n}_chr{c}.trees for the tskit sweep later
keep_all_files: false

end_chr:
  30: {}                # 30 x 100 Mb, constant 1e-8

demographies:
  constant_Ne_10k:
    object: constant_Ne
    path: ibd_sims/demography.py

mating:
  coalescent:
    pedigree_mode: false
```

Check the plan, then create the run directory and per-run `args.yaml`:

```bash
python ibd_sims/experiment.py describe yaml_files/minspan_sweep.yaml
python ibd_sims/experiment.py init     yaml_files/minspan_sweep.yaml
```

This creates `minspan_sweep/constant_Ne_10k__coalescent/` (the `RUN_DIR` below) with its
`args.yaml`.

## Phase 1: simulate and call hap-ibd

```bash
python ibd_sims/experiment.py commands yaml_files/minspan_sweep.yaml     # prints the command
python run.py simulate minspan_sweep/yaml_files/constant_Ne_10k__coalescent.yaml --no-wait
```

Because `iter >= 3`, each iteration is one Slurm job running `sim_workers` chromosomes in parallel.
Memory is `gb × sim_workers` and time is `sim_min × ceil(end_chr / sim_workers)`. For every
iteration and chromosome the job simulates the tree sequence, writes the VCF and `.map`, runs hap-ibd, and
saves `iter{n}_chr{c}.trees` and `.tmrca.pkl`. It then concatenates the chromosomes into
`iter{n}.ibd.gz` and `iter{n}.map` and removes the per-chromosome `.ibd.gz`, `.map` and VCF files.

Check progress, and rerun the same command for anything that failed:

```bash
python ibd_sims/experiment.py status yaml_files/minspan_sweep.yaml
```

Done when 25 `iter{n}.ibd.gz` files and 750 `iter{n}_chr{c}.trees` files exist in `RUN_DIR`.

### hap-ibd timing

Each chromosome task prints `Time to write VCF and call IBD: <seconds> (iter=N, chrom=C)` in its Slurm log under
`RUN_DIR/slurm/`. Average these for the hap-ibd time per chromosome:

```bash
grep -h "Time to write VCF and call IBD" minspan_sweep/constant_Ne_10k__coalescent/slurm/*_log.out \
 | awk '{it=$9; gsub(/[^0-9]/,"",it); s[it]+=$8; n[it]++}
   END {for (i in s) {m=s[i]/n[i]; K++; t+=m; q+=m*m}
        mean=t/K; sd=sqrt((q-K*mean*mean)/(K-1));
        printf "%d iterations, mean %.2f s, SE %.2f s\n", K, mean, sd/sqrt(K)}'
```

Note that 10 chromosomes run in parallel in each job, so this time is measured under that
contention. The tskit time below is measured on a single core, so treat the comparison as
approximate. The tskit time covers only the IBD calling step, and the hap-ibd time also includes
writing the VCF.

## Phase 2: call IBD from the saved trees (tskit, all `min_span`)

One Slurm array task per iteration (`minspan_call.sbatch`):

```bash
#!/bin/bash
#SBATCH --array=1-25
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=minspan_call

python ibd_sims/test/minspan_sweep.py call minspan_sweep/constant_Ne_10k__coalescent \
    --iter $SLURM_ARRAY_TASK_ID
```

```bash
sbatch minspan_call.sbatch
```

For each chromosome, `call` loads `iter{n}_chr{c}.trees`, calls IBD at every `min_span` in
`20000,50000,100000,200000,500000,1000000` (the default; change with `--min-spans`), and caches the
results in `RUN_DIR/sweep/iter{n}_chr{c}_ms{span}.ibd.gz`. It records the wall-clock time of each
call in `RUN_DIR/sweep/iter{n}_chr{c}_timing.tsv`. Cached results are skipped on rerun; use
`--force` to recompute them.

## Phase 3: segment-level comparison

```bash
python ibd_sims/test/minspan_sweep.py summarize minspan_sweep/constant_Ne_10k__coalescent \
    --iters 1-25 --hapibd-sec-per-chr <mean from the awk above> --hapibd-sec-se <SE from the awk above>
```

Writes to `RUN_DIR/sweep/`: `metrics_by_iteration.tsv`, `summary.tsv` (mean and SE across
iterations), and `summary.png`. The printed table shows RMSE and bias of per-pair total IBD (cM),
percent difference in total IBD, segment-count ratios, and seconds per chromosome. Every number is
relative to the 20 kb reference.

## Phase 4: IBDNe

Build one directory per method (hap-ibd and each `min_span`) with genome-level IBD, a shared map,
and an `args.yaml` for the repo's own post-processing, then run IBDNe on each:

```bash
python ibd_sims/test/minspan_sweep.py export minspan_sweep/constant_Ne_10k__coalescent --iters 1-25
bash minspan_sweep/constant_Ne_10k__coalescent/sweep/ne/run_ne.sh
```

`run_ne.sh` runs `python run.py postprocess DIR --no-wait` once per method (IBDNe settings: `mincm 2`,
`trimcm 0.2`, `gmin 1`, `gmax 300`, `nboots 0`, `nits 1000`). Wait for those Slurm jobs to finish, then:

```bash
python ibd_sims/test/minspan_sweep.py summarize-ne minspan_sweep/constant_Ne_10k__coalescent --iters 1-25
```

This reports, per method, the RMSE and bias of log2 Ne against the true Ne (10,000) over
generations 5 to 100, and the RMSE of the curve against the reference method's curve. Outputs are
`ne_metrics_by_iteration.tsv`, `ne_summary.tsv` and `ne_summary.png` in `RUN_DIR/sweep/`.

`export` supports constant-rate runs only, and refuses `end_chr: 22`.

## Phase 5: tables for the paper

```bash
python ibd_sims/test/minspan_sweep.py latex minspan_sweep/constant_Ne_10k__coalescent \
    --which both --out minspan_sweep/tables.tex
```

The tables use `booktabs` and `graphicx` (`\usepackage{booktabs,graphicx}`). Dashes mark quantities
that are not defined (RMSE and bias for the reference, and time for hap-ibd if
`--hapibd-sec-per-chr` was not given).

## Second experiment: `minspan_sweep_100k` (Ne = 100,000)

Same design and the same steps as above, with the demography changed to `constant_Ne100k`
(Ne = 100,000) and everything renamed so the two experiments never share a directory. Replace each
10k command with its 100k counterpart below. Anything not listed here is unchanged.

| | 10k experiment | 100k experiment |
|---|---|---|
| Experiment name | `minspan_sweep` | `minspan_sweep_100k` |
| Experiment yaml | `yaml_files/minspan_sweep.yaml` | `yaml_files/minspan_sweep_100k.yaml` |
| Demography key / object | `constant_Ne_10k` / `constant_Ne` | `constant_Ne_100k` / `constant_Ne100k` |
| Run directory (`RUN_DIR`) | `minspan_sweep/constant_Ne_10k__coalescent` | `minspan_sweep_100k/constant_Ne_100k__coalescent` |
| Call sbatch script | `minspan_call.sbatch` | `minspan_call_100k.sbatch` |
| True Ne | 10,000 | 100,000 |

**Step 0: start clean**

```bash
mv minspan_sweep_100k minspan_sweep_100k_old      # only if an earlier run exists
```

**Step 1: the experiment definition.** `yaml_files/minspan_sweep_100k.yaml` is the 10k file with the
name and the demography changed (all other settings, including 25 iterations, 100 samples and
30 chromosomes of 100 Mb, are identical):

```yaml
experiment: minspan_sweep_100k

# Fixed across all simulations
iter: 25
samples: 100
sim_workers: 10

# Resources (per chromosome; the per-iteration job multiplies these)
gb: 8
sim_min: 30
nthreads: 2

# This experiment
tskit_ibd: false
keep_trees: true
keep_all_files: false

end_chr:
  30: {}

demographies:
  constant_Ne_100k:
    object: constant_Ne100k
    path: ibd_sims/demography.py

mating:
  coalescent:
    pedigree_mode: false
```

```bash
python ibd_sims/experiment.py describe yaml_files/minspan_sweep_100k.yaml
python ibd_sims/experiment.py init     yaml_files/minspan_sweep_100k.yaml
```

**Phase 1: simulate and call hap-ibd**

```bash
python ibd_sims/experiment.py commands yaml_files/minspan_sweep_100k.yaml     # prints the command
python run.py simulate minspan_sweep_100k/yaml_files/constant_Ne_100k__coalescent.yaml --no-wait
python ibd_sims/experiment.py status yaml_files/minspan_sweep_100k.yaml
```

hap-ibd timing:

```bash
grep -h "Time to write VCF and call IBD" minspan_sweep_100k/constant_Ne_100k__coalescent/slurm/*_log.out \
 | awk '{it=$9; gsub(/[^0-9]/,"",it); s[it]+=$8; n[it]++}
   END {for (i in s) {m=s[i]/n[i]; K++; t+=m; q+=m*m}
        mean=t/K; sd=sqrt((q-K*mean*mean)/(K-1));
        printf "%d iterations, mean %.2f s, SE %.2f s\n", K, mean, sd/sqrt(K)}'
```

**Phase 2: call IBD from the saved trees**: `minspan_call_100k.sbatch`

```bash
#!/bin/bash
#SBATCH --array=1-25
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=minspan_call_100k

python ibd_sims/test/minspan_sweep.py call minspan_sweep_100k/constant_Ne_100k__coalescent \
    --iter $SLURM_ARRAY_TASK_ID
```

```bash
sbatch minspan_call_100k.sbatch
```

**Phase 3: segment-level comparison**

```bash
python ibd_sims/test/minspan_sweep.py summarize minspan_sweep_100k/constant_Ne_100k__coalescent \
    --iters 1-25 --hapibd-sec-per-chr <mean from the awk above> --hapibd-sec-se <SE from the awk above>
```

**Phase 4: IBDNe**. The true Ne is different, so pass `--true-ne 100000`:

```bash
python ibd_sims/test/minspan_sweep.py export minspan_sweep_100k/constant_Ne_100k__coalescent --iters 1-25
bash minspan_sweep_100k/constant_Ne_100k__coalescent/sweep/ne/run_ne.sh
# wait for the Slurm jobs, then:
python ibd_sims/test/minspan_sweep.py summarize-ne minspan_sweep_100k/constant_Ne_100k__coalescent \
    --iters 1-25 --true-ne 100000
```

**Phase 5: tables for the paper**

```bash
python ibd_sims/test/minspan_sweep.py latex minspan_sweep_100k/constant_Ne_100k__coalescent \
    --which both --true-ne 100000 --out minspan_sweep_100k/tables.tex
```

The Ne error window stays at generations 5 to 100 (`--gen-min`, `--gen-max`) unless you change it.
With Ne = 100,000 there is less IBD per pair, so expect fewer segments and noisier estimates than
in the 10k experiment.

## Outputs at a glance

| Path (under `RUN_DIR`) | From | Contents |
|---|---|---|
| `iter{n}.ibd.gz`, `iter{n}.map` | Phase 1 | hap-ibd calls, all chromosomes |
| `iter{n}_chr{c}.trees` | Phase 1 | tree sequences |
| `slurm/*_log.out` | Phase 1 | logs, with the hap-ibd timing line |
| `sweep/iter{n}_chr{c}_ms{span}.ibd.gz` | Phase 2 | tskit calls |
| `sweep/iter{n}_chr{c}_timing.tsv` | Phase 2 | tskit call times |
| `sweep/metrics_by_iteration.tsv`, `summary.tsv`, `summary.png` | Phase 3 | segment comparison |
| `sweep/ne/<method>/` | Phase 4 | per-method inputs and IBDNe output |
| `sweep/ne_*.tsv`, `ne_summary.png` | Phase 4 | Ne comparison |

The LaTeX tables are written outside `RUN_DIR`, to `minspan_sweep/tables.tex` (10k) and
`minspan_sweep_100k/tables.tex` (100k).
