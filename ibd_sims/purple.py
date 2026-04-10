import pandas as pd
import itertools as it
import numpy as np
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, List, Union
from multiprocessing import Pool
from functools import partial
import time


def process_dict(pair_dict, ids_to_keep=[]):
    # Get the number of individuals to keep
    n_indv = len(ids_to_keep)

    # True if there are indv to keep
    id_bool = len(ids_to_keep) > 0

    if id_bool:
        mat = np.zeros((300, 300, n_indv+1, 2), dtype=np.int32)
        id_to_idx = {iid: idx for idx, iid in enumerate(ids_to_keep)}
        max_idx = n_indv

        for pair, segments in pair_dict.items():
            if len(segments) == 1:
                continue

            idx1 = id_to_idx.get(pair[0], max_idx)
            idx2 = id_to_idx.get(pair[1], max_idx)

            if len(segments) == 2:
                (i,j,k), (x,y,z) = segments

                k, z = sorted([k, z])

                mat[k, z, [idx1, idx2], 1] += [1, 1]
                mat[k, z, [idx1, idx2], 0] += [i!=x, j!=y]

                continue

            for (i,j,k), (x,y,z) in it.combinations(segments, r=2):

                k, z = sorted([k, z])

                mat[k, z, [idx1, idx2], 1] += [1, 1]
                mat[k, z, [idx1, idx2], 0] += [i!=x, j!=y]

    else:
        mat = np.zeros((300, 300, 2), dtype=np.int32)

        for pair, segments in pair_dict.items():
            if len(segments) == 1:
                continue

            if len(segments) == 2:
                (i,j,k), (x,y,z) = segments

                k, z = sorted([k, z])

                mat[k, z, 1] += 2
                mat[k, z, 0] += sum([i!=x, j!=y])

                continue

            for (i,j,k), (x,y,z) in it.combinations(segments, r=2):

                k, z = sorted([k, z])

                mat[k, z, 1] += 2
                mat[k, z, 0] += sum([i!=x, j!=y])

    return mat


def process_chromosome(filen, ids_to_keep=[], chunksize=1000000):
    pair_dict = {}

    for chunk in pd.read_csv(filen, sep="\\s+", header=None, chunksize=chunksize):

        if len(ids_to_keep) > 0:
            chunk = chunk[(chunk[0].isin(ids_to_keep)) | (chunk[2].isin(ids_to_keep))]

        for row in chunk.itertuples():
            pair = (row[1], row[3])  # Note: indices are offset by 1
            segment = (row[2], row[4], int(row[8]))

            if pair not in pair_dict:
                pair_dict[pair] = [segment]
                continue

            pair_dict[pair].append(segment)


    print("Finished " + filen)

    return pair_dict

def readin_ibd(filen, n_chrom=0, ids_to_keep=[], file_to_write=None, nthreads=8, chunksize=500000):
    t1 = time.time()
    if n_chrom > 0:
        process_with_args = partial(process_chromosome, ids_to_keep=ids_to_keep, chunksize=chunksize)

        with Pool(processes=nthreads) as pool:
            chr_files = [filen.replace("chr1", f"chr{i}") for i in range(1, n_chrom+1)]
            results_list = pool.map(process_with_args, chr_files)

        pair_dict = results_list[0]
        for tmp in results_list[1:]:
            for pair, segments in tmp.items():
                if pair not in pair_dict:
                    pair_dict[pair] = segments
                    continue

                pair_dict[pair].extend(segments) 

    else:
        pair_dict = process_chromosome(filen, ids_to_keep)

    mat = process_dict(pair_dict, ids_to_keep)

    print(time.time()-t1)

    if type(file_to_write) == str:
        np.save(file_to_write, mat)
    
    return mat

# def readin_ibd(filen, n_chrom):
#     df_list = []
#     if "chr" in filen:
#         for c in range(1, n_chrom+1):
#             tmp = pd.read_csv(filen.replace("chr1", f"chr{c}"), sep="\t", header=None)
#             tmp[5] = c
#             tmp[7] = tmp[7].apply(int)
#             df_list.append(tmp)
#         return pd.concat(df_list)
#     df = pd.read_csv(filen, sep="\\s+", header=None)
#     df[7] = df[7].apply(int)
#     return df

def are_purple(seg1, seg2):
    return sum([seg1[0]!=seg1[1], seg2[0]!=seg2[1]])

def find_purple(ibd_df):
    mat = np.zeros((300, 300, 2))
    for _, tmp in ibd_df.groupby([0, 2]):
        for seg1, seg2 in it.combinations(tmp[[1, 3, 7]].sort_values(7).apply(tuple, axis=1).values, r=2):
            mat[seg1[2], seg2[2], 1] += 2 
            mat[seg1[2], seg2[2], 0] += are_purple(seg1, seg2)
    return mat

def purple_matrix(filen, n_chrom):

    ibd = readin_ibd(filen, n_chrom)

    return find_purple(ibd)

def plot(purple_mat, label=None, save_to=None, max_len=20, max_out=False):
    # Calculate the percentage matrix
    percentages = purple_mat[:max_len, :max_len, 0] / purple_mat[:max_len, :max_len, 1]

    if max_out:
        percentages = np.minimum(percentages, 0.5)

    # Create the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(percentages, 
                cmap='YlOrRd',  # You can change the color scheme
                annot=False,     # Show the values in each cell
                fmt='.2%',      # Format as percentage with 2 decimal places
                cbar_kws={'label': 'Purple\nPercentage'})

    if label:
        plt.title(label.replace("\n", " | "))
    else:
        plt.title('Segment Counts Percentage Heatmap')
    plt.xlabel('Length (cM)')
    plt.xlim(2, max_len+1)
    plt.ylim(2, max_len+1)
    plt.ylabel('Length (cM)')

    plt.xticks(rotation=90)
    plt.yticks(rotation=0)

    if save_to:
        plt.savefig(save_to, dpi=500)
        return

    plt.show()

def plot2(
    purple_mat: np.ndarray, 
    label: Optional[str] = None, 
    save_to: Optional[str] = None,
    max_out: bool = False,
    x_labels: Optional[List[str]] = None,
    y_labels: Optional[List[str]] = None
):
    """
    Plot a heatmap of IBD segment comparisons.
    
    Parameters:
    -----------
    purple_mat : np.ndarray
        3D array where purple_mat[:,:,0] contains numerator counts and
        purple_mat[:,:,1] contains denominator counts
    label : str, optional
        Title for the plot
    save_to : str, optional
        Path to save the figure to. If None, the figure is displayed instead
    max_out : bool, default=False
        Whether to cap the maximum proportion at 0.5
    x_labels : List[str], optional
        Custom labels for x-axis ticks. If None, indices are used
    y_labels : List[str], optional
        Custom labels for y-axis ticks. If None, indices are used
    """
    # Get dimensions
    m, n = purple_mat.shape[0], purple_mat.shape[1]
    
    # Calculate the percentage matrix
    percentages = purple_mat[:, :, 0] / np.maximum(purple_mat[:, :, 1], 1)  # Avoid division by zero
    
    # Set zero denominators to NaN to ensure they appear white
    percentages[purple_mat[:, :, 1] == 0] = np.nan
    
    if max_out:
        percentages = np.minimum(percentages, 0.5)
    
    # Create the heatmap
    plt.figure(figsize=(10, 8))
    
    # Create a custom colormap that transitions from white to purple
    cmap = sns.color_palette("Purples", as_cmap=True)
    
    # Generate default labels if needed
    if x_labels is None:
        x_labels = [str(i) for i in range(n)]
    if y_labels is None:
        y_labels = [str(i) for i in range(m)]
    
    # Create the heatmap with mask for NaN values
    # The mask parameter will make NaN values appear white
    mask = np.isnan(percentages)
    
    sns.heatmap(
        percentages, 
        cmap=cmap,      # Purple color scheme
        annot=False,    # Show the values in each cell
        fmt='.2f',      # Format as decimal with 2 decimal places
        xticklabels=x_labels,
        yticklabels=y_labels,
        cbar_kws={'label': 'Proportion'},
        mask=mask       # Apply mask to make NaN values white
    )
    
    if label:
        plt.title(label.replace("\n", " | "))
    else:
        plt.title('IBD Segment Comparison Heatmap')
    
    plt.xlabel('Segment Length j')
    plt.ylabel('Segment Length i')
    
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    
    if save_to:
        plt.savefig(save_to, dpi=500, bbox_inches='tight')
    else:
        plt.show()

if __name__ == "__main__":

    filen = sys.argv[1]

    readin_ibd(filen, n_chrom=int(sys.argv[2]), file_to_write=sys.argv[3])

    # mat = purple_matrix(filen, int(sys.argv[2]) if len(sys.argv)>2 else 30)

    # if len(sys.argv) > 2:

    #     np.save(sys.argv[3], mat)
    