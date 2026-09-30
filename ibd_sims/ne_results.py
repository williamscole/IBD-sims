"""
ne_results.py -- load and summarise saved Ne estimates (Ne_data.pkl).

Deliberately light: depends only on numpy, pandas and the standard library, so it
can be imported from a figures/analysis environment that does not have msprime,
submitit, etc.  Do NOT add imports of other ibd_sims modules here.
"""

import pickle
from pathlib import Path

import numpy as np
import pandas as pd


def load_ne_data(exp_dir, pickle_name="plots/Ne_data.pkl", halve_hapne=True):
    """
    Load [exp_dir]/plots/Ne_data.pkl (written by plot_Ne.py --save-pickle).

    HapNe reports Ne in haploid units; IBDNe and the simulated truth are diploid.
    With halve_hapne=True (default) every HapNe record's NE column is divided by 2.
    The payload gets hapne_halved=True afterwards, so loading a pickle regenerated
    with the fix applied at the source cannot halve twice.
    """
    with open(Path(exp_dir) / pickle_name, "rb") as fh:
        data = pickle.load(fh)

    if halve_hapne and not data.get("hapne_halved", False):
        for demo_payload in data["demos"].values():
            for rec in demo_payload["records"]:
                if "hapne" not in rec["method"]:
                    continue
                for df in rec["dfs"]:
                    assert isinstance(df, pd.DataFrame) and {"GEN", "NE"} <= set(df.columns)
                    df["NE"] = df["NE"] / 2
        data["hapne_halved"] = True

    return data
