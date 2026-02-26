import numpy as np
import pandas as pd
import sys

if __name__ == "__main__":

    segments = pd.read_csv(sys.argv[-2], delim_whitespace=True, header=None).head(10000)

    samples = int(sys.argv[-1])

    df_samples = list(set(segments[0]) | set(segments[2]))

    keep_samples = np.random.choice(df_samples, samples, replace=False)

    segments = segments[segments.apply(lambda x: x[0] in keep_samples and x[2] in keep_samples, axis=1)]

    import pdb; pdb.set_trace()