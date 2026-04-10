import os
import re
from pathlib import Path
from datetime import datetime
import sys

def checked_completed(path, array_n):
    """
    Finds the newest Slurm output file matching the pattern:
    f"{path}/slurm/simulation_array{array_n}_[any jobid].out"
    
    Parameters:
    -----------
    array_n : int or str
        The array number to search for
    path : str
        Base path where the slurm directory is located
        
    Returns:
    --------
    str or None
        Path to the newest matching file, or None if no matches found
    """
    slurm_dir = Path(f"{path}/slurm")
    
    
    # Pattern to match: simulation_array{array_n}_[any jobid].out
    pattern = re.compile(f"simulation_array{array_n}_\\d+\\.out$")
    
    matching_files = []
    
    # Find all matching files
    for file_path in slurm_dir.glob(f"simulation_array{array_n}_*.out"):
        if pattern.match(file_path.name):
            # Get modification time and add to list
            mod_time = file_path.stat().st_mtime
            matching_files.append((file_path, mod_time))
    
    # Return None if no matching files
    if not matching_files:
        return False
    
    # Sort by modification time (newest last) and get the newest file
    matching_files.sort(key=lambda x: x[1])
    newest_file = matching_files[-1][0]

    try:
        i = open(newest_file).readlines()[-1].strip()
    except:
        return False


    return i == "Success!"

def get_arrays_to_run(path, n_arrays, max_jobs=800):

    to_add = [i for i in range(1, n_arrays+1) if not checked_completed(path, i)]

    if len(to_add) == 0:
        return f"1-{n_arrays}"

    units = []
    unit = [to_add[0]]
    for index in range(1, len(to_add)):
        # Not consecutive
        if to_add[index] != (unit[-1] + 1):
            units.append(unit)
            unit = [to_add[index]]
        else:
            unit.append(to_add[index])

    units.append(unit)

    runs = [[]]
    for unit in units:
        cur = sum([len(i) for i in runs[-1]]) # Number of jobs in the current array
        if cur < max_jobs: # If number of jobs is less than the max number of jobs, may be able to add all
            can_add = max_jobs - cur # Determine how many we can add
            if len(unit) <= can_add: # We can add all of them
                runs[-1].append(unit)
            else: # We can only add some of them
                runs[-1].append(unit[:can_add]) # Add the rest
                unit = unit[can_add:]
                while len(unit) > 0:
                    runs.append([unit[:max_jobs]])
                    unit = unit[max_jobs:]
        else:
            if len(unit) <= max_jobs:
                runs.append([unit])
            else:
                while len(unit) > 0:
                    runs.append([unit[:max_jobs]])
                    unit = unit[max_jobs:]

    # import pdb; pdb.set_trace()

    for run in runs:
        print(",".join([str(i[0]) if len(i)==1 else f"{i[0]}-{i[-1]}" for i in run]))





if __name__ == "__main__":


    get_arrays_to_run(sys.argv[1], int(sys.argv[2]))
