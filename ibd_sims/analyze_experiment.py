import re
import pandas as pd
from pathlib import Path
import sys

def _load_ibdne_file(ne_file: Path) -> pd.DataFrame:
    """Parse an IBDNe .ne output file into a DataFrame.

    Expected format (tab-separated):
        GEN     NE      LWR-95%CI       UPR-95%CI
        0       3.36E4  3.36E4          3.36E4
    """
    df = pd.read_csv(ne_file, sep="\t")
    df.columns = ["gen", "ne", "lwr_95ci", "upr_95ci"]
    return df


def _load_hapne_file(hapne_csv: Path) -> pd.DataFrame:
    """Parse a HapNe hapne.csv output file into a DataFrame.

    Expected format (comma-separated):
        TIME,Q0.025,Q0.25,Q0.5,Q0.75,Q0.975
        1,13457,14672,14766,16583,19501
    """
    df = pd.read_csv(hapne_csv)
    df.columns = ["time", "ne_q025", "ne_q25", "ne_q50", "ne_q75", "ne_q975"]
    return df


def _load_postprocess_args(exp_dir: Path) -> dict:
    """
    Load postprocess.tsv and return a per-method lookup DataFrame mapping
    rep (e.g. "001") to its postprocessing arg columns.

    Parameters
    ----------
    exp_dir : Path
        Experiment directory containing postprocess.tsv.

    Returns
    -------
    dict keyed by method name (e.g. "ibdne", "hapne_ibd"), each value a
    DataFrame with columns [rep, <arg1>, <arg2>, ...].
    Returns an empty dict if postprocess.tsv is not found.
    """
    tsv_path = exp_dir / "postprocess.tsv"
    if not tsv_path.exists():
        return {}

    pp = pd.read_csv(tsv_path, sep="\t")

    # directory column is like "ibdne/001" or "hapne_ibd/002"
    pp["rep"] = pp["directory"].str.split("/").str[-1]

    method_args = {}
    for method, group in pp.groupby("name"):
        # Collect the union of arg column names declared for this method
        arg_cols = set()
        for args_str in group["args"].dropna():
            arg_cols.update(a.strip() for a in args_str.split(","))
        arg_cols = sorted(arg_cols)

        keep_cols = ["rep"] + [c for c in arg_cols if c in group.columns]
        method_args[method] = group[keep_cols].reset_index(drop=True)

    return method_args


def load_experiment_results(exp_dir):
    """
    Load Ne estimates from ibdne and hapne_ibd postprocessing runs for an experiment.

    Also loads postprocess.tsv (if present) to annotate each row with the
    postprocessing args that correspond to that rep directory (e.g. filter,
    filtersamples).

    Parameters
    ----------
    exp_dir : str or Path
        Path to the experiment directory. Must contain yaml_files/yaml_files.txt
        and subdirectories named by yaml prefix (e.g. constant_Ne_10k__DTWF_di).

    Returns
    -------
    dict with keys "ibdne" and "hapne_ibd", each a pandas DataFrame:

        ibdne columns:
            demo, rep, iter, <postprocess args...>, gen, ne, lwr_95ci, upr_95ci

        hapne_ibd columns:
            demo, rep, iter, <postprocess args...>, time, ne_q025, ne_q25, ne_q50, ne_q75, ne_q975
    """
    exp_dir = Path(exp_dir)

    yaml_files_txt = exp_dir / "yaml_files" / "yaml_files.txt"
    with open(yaml_files_txt) as f:
        demos = [line.strip().replace(".yaml", "") for line in f if line.strip()]

    ibdne_rows = []
    hapne_rows = []

    for demo in demos:
        demo_dir = exp_dir / demo

        # --- ibdne ---
        # Structure: [demo]/ibdne/001/iter5.ne
        ibdne_dir = demo_dir / "ibdne"
        if ibdne_dir.exists():
            for rep_dir in sorted(ibdne_dir.iterdir()):
                if not rep_dir.is_dir() or not re.fullmatch(r"\d{3}", rep_dir.name):
                    continue
                rep = rep_dir.name
                for ne_file in sorted(rep_dir.glob("iter*.ne")):
                    m = re.search(r"iter(\d+)", ne_file.stem)
                    if not m:
                        continue
                    iter_num = int(m.group(1))
                    try:
                        df = _load_ibdne_file(ne_file)
                    except Exception as e:
                        print(f"  [skip] could not load {ne_file}: {e}")
                        continue
                    df.insert(0, "demo", demo)
                    df.insert(1, "rep", rep)
                    df.insert(2, "iter", iter_num)
                    ibdne_rows.append(df)

        # --- hapne_ibd ---
        # Structure: [demo]/hapne_ibd/001/iter1/HapNe/hapne.csv
        hapne_dir = demo_dir / "hapne_ibd"
        if hapne_dir.exists():
            for rep_dir in sorted(hapne_dir.iterdir()):
                if not rep_dir.is_dir() or not re.fullmatch(r"\d{3}", rep_dir.name):
                    continue
                rep = rep_dir.name
                for iter_dir in sorted(rep_dir.iterdir()):
                    if not iter_dir.is_dir():
                        continue
                    m = re.fullmatch(r"iter(\d+)", iter_dir.name)
                    if not m:
                        continue
                    iter_num = int(m.group(1))
                    hapne_csv = iter_dir / "HapNe" / "hapne.csv"
                    if not hapne_csv.exists():
                        print(f"  [skip] not found: {hapne_csv}")
                        continue
                    try:
                        df = _load_hapne_file(hapne_csv)
                    except Exception as e:
                        print(f"  [skip] could not load {hapne_csv}: {e}")
                        continue
                    df.insert(0, "demo", demo)
                    df.insert(1, "rep", rep)
                    df.insert(2, "iter", iter_num)
                    hapne_rows.append(df)

    ibdne_df = (
        pd.concat(ibdne_rows, ignore_index=True) if ibdne_rows else pd.DataFrame()
    )
    hapne_df = (
        pd.concat(hapne_rows, ignore_index=True) if hapne_rows else pd.DataFrame()
    )

    # Annotate with postprocessing args from postprocess.tsv
    pp_args = _load_postprocess_args(exp_dir)

    if not ibdne_df.empty and "ibdne" in pp_args:
        ibdne_df = ibdne_df.merge(pp_args["ibdne"], on="rep", how="left")

    if not hapne_df.empty and "hapne_ibd" in pp_args:
        hapne_df = hapne_df.merge(pp_args["hapne_ibd"], on="rep", how="left")

    return {"ibdne": ibdne_df, "hapne_ibd": hapne_df}


def main():

    exp_dir = sys.argv[1]

    dfs = load_experiment_results(exp_dir)

    import pdb; pdb.set_trace()

if __name__ == "__main__":
    main()
