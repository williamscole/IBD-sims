import pickle as pkl
import yaml
import sys
import os
import numpy as np
import pandas as pd
import ast

from simulations import TMRCA

def fix_zero(proportion, arr):
    out = []
    for i,j in zip(proportion, arr):
        if round(i, 3) > 0:
            out.append(str(j))
    return out

def load_tmrca(file_path, sep='\t'):
    """
    Load chromosome data from a file with list-like strings and convert them to appropriate types.
    
    Parameters:
    -----------
    file_path : str
        Path to the data file
    sep : str, default='\t'
        Separator used in the file
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with columns converted to appropriate types
    """
    # Read the file with correct separator
    df = pd.read_csv(file_path, sep=sep)

    if type(df.iloc[0]["chromosome"]) == str:
        df = df[df.chromosome!="chromosome"]
        df["chromosome"] = df.chromosome.apply(int)
        
    # Function to safely convert string representations of lists to actual lists
    def safe_eval(x):
        if pd.isna(x):
            return []
        if isinstance(x, str):
            if x.startswith('[') and x.endswith(']'):
                try:
                    return ast.literal_eval(x)
                except (ValueError, SyntaxError):
                    return [x]  # Return as single item list if evaluation fails
            else:
                # Handle the case where the value is a string but not a list representation
                return [x]
        return x
    
    # Process each column based on its content
    for col in df.columns:
        if col == 'chrom_index':
            # Convert to integer
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
        elif col == 'chromosome':
            # Keep as is or convert if needed
            pass
        else:
            # For columns that contain list representations
            df[col] = df[col].apply(safe_eval)
            
            # Further processing: convert list elements to appropriate types
            def process_list_elements(lst):
                if not lst or not isinstance(lst, list):
                    return lst
                
                # Try to convert string numbers in lists to float or int
                processed = []
                for item in lst:
                    if isinstance(item, str) and item.replace('.', '', 1).isdigit():
                        # Check if it's a float or int
                        if '.' in item:
                            processed.append(float(item))
                        else:
                            processed.append(int(item))
                    else:
                        processed.append(item)
                return processed
            
            df[col] = df[col].apply(process_list_elements)
    
    return df


# Example usage:
# df = load_chromosome_data('chromosome_data.txt')
# print(df)

def concat_tmrca(path, iter_n, end_chr):

    files = [f"{path}/iter{iter_n}_chr{chrom}.tmrca.pkl" for chrom in range(1, end_chr + 1)]

    assert np.mean([os.path.exists(f) for f in files]) == 1

    out_df = []

    file_to_write =  f"{path}/iter{iter_n}.tmrca.gz"

    for chrom, f in enumerate(files):
        with open(f, "rb") as pklf:
            tmrca = pkl.load(pklf)

            max_index = max(tmrca.tmrca.keys())

            assert len(tmrca.tmrca.keys()) == (max_index + 1)

            df = pd.DataFrame({"chrom_index": np.arange(0, max_index+1)})
            df["chromosome"] = chrom + 1

            df["proportion"] = df["chrom_index"].apply(lambda x: tmrca.get_info(x, "proportion"))

            df["tmrca"] = df.apply(lambda x: fix_zero(x.proportion, tmrca.get_info(x.chrom_index, "tmrca")), axis=1)

            mrcas = df.apply(lambda x: fix_zero(x.proportion, tmrca.get_info(x.chrom_index, "mrca")), axis=1)
            indv = mrcas.apply(lambda x: [tmrca.node_to_indv[int(i)] for i in x])
            indv_name = indv.apply(lambda x: [tmrca.indv_to_name[i] for i in x])

            df["mrca"] = indv_name
            df["proportion"] = df.apply(lambda x: fix_zero(x.proportion, tmrca.get_info(x.chrom_index, "proportion")), axis=1)

            if chrom == 0:
                df.to_csv(file_to_write, sep="\t", index=False, header=True, compression='gzip', mode='w')
            else:
                df.to_csv(file_to_write, sep="\t", index=False, header=False, compression='gzip', mode='a')

    print("Wrote: " + file_to_write)


if __name__ == "__main__":

    path = sys.argv[1]
    iter_n = sys.argv[2]

    end_chr = yaml.safe_load(open(f"{path}/args.yaml"))["end_chr"]

    concat_tmrca(path, iter_n, end_chr)