import yaml
import itertools as it
import sys
import os
import pandas as pd
import pickle as pkl

def load_default():
    return yaml.safe_load(open("yaml_files/args.yaml"))

def write_yaml(dir_name, params, n=1):

    yargs = load_default()

    yargs["dir_name"] = dir_name

    param_names = list(params.keys())
    param_vals = params.values()

    for combo in it.product(*param_vals):
        out_yargs = yargs.copy()
        for index, val in enumerate(combo):
            if type(val) == dict:
                for key, val in val.items():
                    out_yargs[param_names[index]][key] = val
            else:
                out_yargs[param_names[index]] = val

        out_yargs = get_label(out_yargs)

        out_yargs = rules(out_yargs)

        with open(f"yaml_files/arg{n}.yaml", "w") as yout:
            yaml.dump(out_yargs, yout, default_flow_style=False)
            print(f"yaml_files/arg{n}.yaml")
            n += 1

    return n

def update_dict(d, path_value):
    # The last element in path_value is the value to set
    path = path_value[:-1]
    value = path_value[-1]
    
    # Start with the main dictionary
    current = d
    
    # Traverse through all but the last key
    for i, key in enumerate(path[:-1]):
        # If the key doesn't exist or its value isn't a dictionary,
        # create a new dictionary
        if key not in current or not isinstance(current[key], dict):
            current[key] = {}
        current = current[key]
    
    # Set the final value
    current[path[-1]] = value
    
    return d

def get_label(yargs):

    n = yargs["samples"]

    custom_demo = yargs["custom_demo"]["object"]
    custom_sim = yargs["custom_sim"]["object"]

    filtered = yargs["filtersamples"]
    ibd_filter = yargs["ibd_filter"]

    label = []

    if custom_demo:
        label.append(custom_demo)
    elif custom_sim:
        label.append(custom_sim)
    
    label.append(f"n={n}")

    if yargs["pedigree"]["pedigree_mode"]:
        mating = yargs["pedigree"]["mating"]
        label.append(f"DTWF ({mating})")
    elif type(yargs["custom_sim"]["object"]) == str:
        label.append(yargs["custom_sim"]["object"])
    else:
        label.append("coalescent")

    if ibd_filter:
        label.append(ibd_filter)
    else:
        label.append("filtered" if yargs["filtersamples"] else "unfiltered")

    yargs["label"] = "\n".join(label)

    return yargs


def change_in_place(paths, changes):
    def to_bool(yarg):
        if yarg == "false":
            return False
        if yarg == "true":
            return True
        if yarg == "null":
            return None
        try:
            return int(yarg)
        except:
            return yarg

    def get_changes(changes):
        tmp = []
        for change in changes:
            change = change.split(":")
            change[-1] = to_bool(change[-1])
            tmp.append(change)
        return tmp

    changes = get_changes(changes)

    for path in paths:

        cur_yaml = yaml.safe_load(open(f"{path}/args.yaml"))

        for change in changes:
            cur_yaml = update_dict(cur_yaml, change)

        cur_yaml = rules(cur_yaml)

        cur_yaml = get_label(cur_yaml)

        with open(f"{path}/args.yaml", "w") as yout:
            yaml.dump(cur_yaml, yout, default_flow_style=False)
            print(f"Wrote: {path}/args.yaml")


def rules(args):
    if args["simulate_only"]:
        args["run_ibdne"] = False
        args["purple_nodes"] = False

    if args["simulate_only"] or args["run_ibdne"] == False:
        args["nthreads"] = 1

    if args["post_process_only"] and args["run_ibdne"] == False and args["purple_nodes"]:
        args["gb"] = 8

    if args["custom_demo"]["object"] == "constant_Ne100k" or args["custom_demo"]["object"] == "ooa2":
        if args["samples"] == 10000:
            args["gb"] = 8
            args["sim_time"] = "--time=6:00:00"
            args["iter"] = 25
            args["end_chr"] = 30
            args["nthreads"] = 8
            args["ibdne_time"] = "--time=3:00:00"
        if args["samples"] == 1000:
            args["sim_time"] = "--time=1:30:00"

    if args["ibd_filter"]:
        args["filtersamples"] = False

    if args["nboots"] > 80:
        cur_hour = args["ibdne_time"].split("=")[1][0]
        args["ibne_time"] = args["ibdne_time"].replace(f"={cur_hour}", f"={int(cur_hour)+2}")

    return args


if __name__ == "__main__":

    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):

        paths = [i for i in sys.argv[1:] if os.path.exists(i)]

        changes = [i for i in sys.argv[1:] if not os.path.exists(i)]

        change_in_place(paths, changes)

        sys.exit()

    n = write_yaml(sys.argv[1],
    {
        "custom_demo": [{"path": "demography.py", "object": "constant_Ne"},
                        {"path": "demography.py", "object": "constant_Ne100k"}],
        "pedigree": [{"pedigree_mode": True, "mating": "di"},
                     {"pedigree_mode": True, "mating": "mono"}],
        "samples": [1000],
        "iter": [50],
        "nboots": [80],
        "filtersamples": [False],
        "run_ibdne": [True],
        "purple_nodes": [True],
        "simulate_only": [False],
        "post_process_only": [False]
    }
)

    n = write_yaml(sys.argv[1],
    {
        "custom_demo": [{"path": "demography.py", "object": "ooa2"}],
        "pedigree": [{"pedigree_mode": True, "mating": "di"},
                     {"pedigree_mode": True, "mating": "mono"}],
        "samples": [2000],
        "iter": [50],
        "nboots": [80],
        "filtersamples": [False],
        "run_ibdne": [True],
        "purple_nodes": [True],
        "simulate_only": [False],
        "post_process_only": [False]
    },
    n
)


    n = write_yaml(sys.argv[1],
        {
            "custom_sim": [{"path": "load_quebec.py", "object": "load_random_10000"}],
            "ibdne_time": ["--time=1:00:00"],
            "sim_time": ["--time=30:00"],
            "samples": [10000],
            "nthreads": [8],
            "nboots": [80],
            "gb": [32],
            "end_chr": [22],
            "iter": [1],
            "filtersamples": [False],
            "run_ibdne": [True],
            "purple_nodes": [True],
            "simulate_only": [False],
            "post_process_only": [False],
        },
        n
    )

    n = write_yaml(sys.argv[1],
    {
        "custom_demo": [{"path": "demography.py", "object": "ashkenazi"}],
        "pedigree": [{"pedigree_mode": True, "mating": "di"}],
        "samples": [1000],
        "iter": [50],
        "nboots": [80],
        ""
        "filtersamples": [False],
        "run_ibdne": [True],
        "purple_nodes": [True],
        "simulate_only": [False],
        "post_process_only": [False]
    },
    n
)


#     n = write_yaml(sys.argv[1],
#     {
#         "custom_demo": [{"path": "demography.py", "object": "constant_Ne100k"}],
#         "pedigree": [{"pedigree_mode": True, "mating": "di"}],
#         "samples": [10000],
#         "iter": [25],
#         "filtersamples": [False],
#         "run_ibdne": [True],
#         "purple_nodes": [True],
#         "simulate_only": [False],
#         "post_process_only": [False]
#     },
#     n
# )

#     n = write_yaml(sys.argv[1],
#     {
#         "custom_demo": [{"path": "demography.py", "object": "ooa2"}],
#         "pedigree": [{"pedigree_mode": True, "mating": "di"}],
#         "samples": [1000, 10000],
#         "iter": [1],
#         "filtersamples": [False],
#         "run_ibdne": [True],
#         "purple_nodes": [True],
#         "simulate_only": [False],
#         "post_process_only": [False]
#     },
#     n
# )

    
#     write_yaml(sys.argv[1],
#     {
#         "custom_demo": [{"path": "demography.py", "object": "constant_Ne100k"}],
#         "pedigree": [{"pedigree_mode": True, "mating": "di"}, {"pedigree_mode": False}],
#         "samples": [1000],
#         "iter": [10],
#         "filtersamples": [False],
#         "run_ibdne": [True],
#         "purple_nodes": [True],
#         "simulate_only": [False],
#         "post_process_only": [False],
#         "end_chr": [30],
#         "nthreads": [8],
#         "sim_time": ["--time=3:00:00"]
#     }
# )
#     write_yaml(sys.argv[1],
#     {
#         "custom_demo": [{"path": "demography.py", "object": "constant_Ne"},
#                         {"path": "demography.py", "object": "himba"},
#                         {"path": "demography.py", "object": "ooa2"}],
#         "pedigree": [{"pedigree_mode": True, "mating": "di"}],
#         "samples": [1000],
#         "iter": [10],
#         "filtersamples": [False],
#         "run_ibdne": [False],
#         "purple_nodes": [True],
#         "simulate_only": [False],
#         "post_process_only": [False]
#     }
# )