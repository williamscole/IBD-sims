import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import sys
import yaml
import os

from simulations import DemographicSetup

def plot_ne_estimates(data_dict, truth_df=None, figsize=(12, 8), xlim=(0, 100), 
                     log_scale=True, colors=None):
    """
    Plot effective population size estimates with confidence intervals.
    
    Parameters:
    -----------
    data_dict : dict
        Dictionary where keys are labels and values are lists of DataFrames.
        Each DataFrame should have 'GEN' and 'NE' columns.
    truth_df : pd.DataFrame, optional
        DataFrame with 'GEN' and 'NE' columns representing the true values
    figsize : tuple, optional
        Figure size (width, height)
    xlim : tuple, optional
        x-axis limits (min, max)
    log_scale : bool, optional
        Whether to use log scale for y-axis
    colors : list, optional
        List of colors for each dataset. If None, uses a predefined color palette
        
    Returns:
    --------
    fig : matplotlib.figure.Figure
        The figure object
    ax : matplotlib.axes.Axes
        The axes object
    """
    # Set up the plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Define marker styles for different lines
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', 'h', '8']
    # Define line styles
    line_styles = ['-', '--', '-.', ':']
    
    # Generate colors if not provided
    if colors is None:
        # Create a color palette with distinct, visually pleasing colors
        n_colors = len(data_dict)
        colors = sns.color_palette("deep", n_colors)  # Changed to a more distinct palette

    title_strs = []
    tmp = []
    for label in data_dict.keys():
        tmp.extend(label.split("\n"))
    for i in tmp:
        if tmp.count(i) == len(data_dict.keys()):
            title_strs.append(i)
    plot_title = " | ".join(list(set(title_strs)))

    
    # Plot each dataset
    for i, ((label, dfs), color) in enumerate(zip(data_dict.items(), colors)):
        # Convert list of DataFrames to 3D array
        data_array = np.array([df['NE'].values for df in dfs])

        # Calculate statistics
        mean = np.mean(data_array, axis=0)
        
        # Get generations (using first DataFrame as reference)
        generations = dfs[0]['GEN']
        
        # Calculate marker positions - stagger them between lines
        n_markers = 15  # Number of markers per line
        marker_indices = np.linspace(0, len(generations)-1, n_markers, dtype=int)
        # Offset marker positions for each line to avoid overlap
        offset = i / len(data_dict) * (len(generations) // n_markers)
        marker_indices = [int(idx + offset) % len(generations) for idx in marker_indices]
        
        # Create marker array
        markevery = [False] * len(generations)
        for idx in marker_indices:
            if idx < len(markevery):
                markevery[idx] = True
        
        # Plot mean with markers and line
        marker = markers[i % len(markers)]  # Cycle through markers
        linestyle = line_styles[i % len(line_styles)]  # Cycle through line styles
        
        ax.plot(generations, mean, 
               label="\n".join([i for i in label.split("\n") if i not in title_strs]), 
               color=color, 
               linewidth=2.5,           # Increased line width
               linestyle=linestyle,)     # Different line style for each line

        if len(data_array) > 1:
            percentile_5 = np.percentile(data_array, 5, axis=0)
            percentile_95 = np.percentile(data_array, 95, axis=0)
            ax.fill_between(generations, percentile_5, percentile_95,
                          color=color, alpha=0.15)  # Reduced alpha for less visual noise
    
    # Plot truth data if provided
    if truth_df is not None:
        ax.plot(truth_df['GEN'], truth_df['NE'], 
                color='black',
                linestyle=(0, (3, 5, 1, 5, 1, 5)),          # Solid line for truth
                label='Truth',
                linewidth=2)
                # marker='x',
                # markevery=max(1, len(truth_df)//15),
                # markersize=10,
                # markeredgewidth=2)
                
        if np.std(truth_df["NE"]) == 0:
            ax.axvline(x=np.log(truth_df.iloc[0]["NE"]) / np.log(2), 
                      linestyle="--", color="black", linewidth=1.5,
                      label=r"$log_2N_e$")
            ax.axvline(x=1.77 * np.log(truth_df.iloc[0]["NE"]) / np.log(2), 
                      linestyle=":", color="black", linewidth=2,
                      label=r"$1.77 \times log_2N_e$")

    # ax.set_xlim(xlim)
    ax.set_xlim(0, 50)
    # Customize the plot
    if log_scale:
        ax.set_yscale('log')
    
    ax.set_title(plot_title,
                pad=20)
    ax.set_xlabel('Generation')
    ax.set_ylabel('Effective Population Size')
    ax.grid(True, alpha=0.2)  # Reduced grid opacity
    
    # Improve legend
    ax.legend(bbox_to_anchor=(1.05, 1), 
             loc='upper left',
             borderaxespad=0.,
             framealpha=1.0)  # Solid legend background
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    return fig, ax

def get_truth(args):

    demo = DemographicSetup.create(args)

    gen_arr = np.arange(0, args["gmax"]+1)

    ne_arr = demo.debug().population_size_trajectory(gen_arr)

    return pd.DataFrame({"GEN": gen_arr, "NE": ne_arr[:,0]})

def meta_analysis(iter1_path, n_iter, skip_error=True):

    mat_arr = []
    for i in range(1, n_iter+1):
        if not os.path.exists(iter1_path.replace("iter1", f"iter{i}")):
            print("IBDNe file %s not found." % iter1_path.replace("iter1", f"iter{i}"))
            continue
        df = pd.read_csv(iter1_path.replace("iter1", f"iter{i}"), sep="\t")
        df["NE"] = df["NE"].apply(lambda x: x + np.random.normal(x, 0.01*x))
        df["LWR-95%CI"] = df["LWR-95%CI"].apply(lambda x: np.random.normal(x, (0.01 if i != 5 else 5)*x))
        df["UPR-95%CI"] = df["UPR-95%CI"].apply(lambda x: np.random.normal(x, (0.01 if i != 5 else 5)*x))
        mat = df[["NE", "LWR-95%CI", "UPR-95%CI"]].values
        mat_arr.append(mat)

    mat = np.stack(mat_arr, axis=0)

    se = ((mat[:,:,2] - mat[:,:,1]) / (2*1.96))**2 # Get SE of each iteration at each estimate 
    between_var = np.tile(np.std(mat[:,:,0], axis=0)**2, (se.shape[0], 1)) # Variance of NE estimates across iterations
    weights = (1/(se + between_var)) # Weighting by the SE of each iteration and the variance across iterations
    pooled_se = np.sqrt(1 / np.sum(weights, axis=0)) # Compute the pooled SE

    pooled_mean = np.sum(mat[:,:,0]*weights, axis=0) / np.sum(weights, axis=0)
    pooled_lwr = pooled_mean - (pooled_se*1.96)
    pooled_upr = pooled_mean + (pooled_se*1.96)

    pooled_mat = np.stack((pooled_mean, pooled_lwr, pooled_upr), axis=1)

    return pd.DataFrame(pooled_mat, columns=["NE", "LWR-95%CI", "UPR-95%CI"])
    

def load_dfs(path):
    yargs = yaml.safe_load(open(f"{path}/args.yaml"))

    df_list = []
    for i in range(1, yargs["iter"]+1):
        if os.path.exists(f"{path}/iter{i}.ne"):
            df_list.append(pd.read_csv(f"{path}/iter{i}.ne", sep="\t"))
        else:
            print(f"Missing iteration {i} for arg file: {path}")

    return yargs, yargs["label"], df_list


def find_ibdne_directories(directory_path):
    ibdne_directories = []
    
    for item in os.listdir(directory_path):
        item_path = os.path.join(directory_path, item)
        
        if os.path.isdir(item_path) and item.endswith("_ibdne"):
            ibdne_directories.append(directory_path + "/" + item)
    
    return ibdne_directories

def plot(path, exclude=[]):

    data_dict = {}

    ibdne_dir = find_ibdne_directories(path)

    for i in exclude:
        if i in ibdne_dir:
            ibdne_dir.remove(i)

    for i in ibdne_dir:

        yargs, label, dfs = load_dfs(i)

        if len(dfs) > 0:
            data_dict[label] = dfs
        else:
            print(f"No data found for path {i}")

    if yargs["custom_sim"]["path"]:
        truth_df = None

    else:
        truth_df = get_truth(yargs)

    fig, ax = plot_ne_estimates(data_dict, truth_df=truth_df)

    outputf = f"{path}/plot.png"

    plt.savefig(outputf, dpi=600)

    print(f"Saved file: {outputf}")


# Example usage:
if __name__ == "__main__":

    plot(sys.argv[1], sys.argv[2:] if len(sys.argv)>2 else [])



    # # Example of how to prepare the data
    # Ne = 10000
    
    # def load_data(path_pattern, n_iterations=50):
    #     return [pd.read_csv(f"{path_pattern.format(i=i, Ne=Ne)}", 
    #                        sep="\\s+")
    #             for i in range(1, n_iterations + 1)]
    
    # # Load datasets
    # data_dict = {
    #     "DTWF": load_data("pedigree/iter{i}_Ne{Ne}_n1000_pedigree.ne"),
    #     "Hudson": load_data("coalescent/iter{i}_Ne{Ne}_n1000_coalescent.ne"),
    #     "Monog DTWF": load_data("monogamous/iter{i}_Ne{Ne}_n1000_pedigree.ne")
    # }
    
    # # Create truth DataFrame
    # truth_df = pd.DataFrame({
    #     'GEN': data_dict['DTWF'][0]['GEN'],
    #     'NE': Ne
    # })
    
    # # Create the plot
    # fig, ax = plot_ne_estimates(data_dict, truth_df=truth_df)
    # plt.show()

