import pandas as pd
import sys
import time
import polars as pl

if __name__ == "__main__":

    df = pd.read_csv(sys.argv[-1], sep="\\s+", header=None)

    counts = {}
    for _, pair_df in df.groupby([0, 2]):
        l = int(pair_df[7].sum())
        counts[l] = counts.get(l, 0) + 1

    print("l", "count")
    for i,j in counts.items():
        print(f"{i} {j}")
