import pandas as pd
import sys
import numpy as np

DEGREES = {
    "11th": (2**-10.5, 2**-9.9),
    "10th": (2**-9.5, 2**-8.5),
    "9th": (2**-8.5, 2**-7.5),
    "8th": (2**-7.5, 2**-6.5),
    "7th": (2**-6.5, 2**-5.5),
    "6th": (2**-5.5, 2**-4.5),
    "5th": (2**-4.5, 2**-3.5),
    "4th": (2**-3.5, 2**-2.5),
    "3rd": (2**-2.5, 2**-1.5),
    "2nd": (2**-1.5, 2**-.5),
    "1st": (2**-.5, np.inf)
}

if __name__ == "__main__":

    df = pd.read_csv(sys.argv[-1], delim_whitespace=True, header=None)

    samples = list(set(df[0].values) | set(df[2].values))
    degrees = {d: index for index, d in enumerate(DEGREES.keys())}

    mat = np.zeros((len(samples), len(degrees)))

    for (id1, id2), pair_df in df.groupby([0, 2]):

        k = pair_df[7].sum() / 3000

        for d, (i, j) in DEGREES.items():
            if i < k <= j:
                idx1 = samples.index(id1)
                idx2 = samples.index(id2)
                dx = degrees[d]
                mat[idx1, dx] += 1
                mat[idx2, dx] += 1

    np.save(sys.argv[-1].replace(".ibd", ".mat.npy"), mat)
    print(f"Saved to {sys.argv[-1].replace('.ibd', '.mat.npy')}")



