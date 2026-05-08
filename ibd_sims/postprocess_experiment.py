import itertools as it
import pandas as pd
import yaml
from pathlib import Path

def load_yaml(yaml_file):
    with open(yaml_file, "r") as yamlf:
        exp_args = yaml.safe_load(yamlf)
    return exp_args

def process_yaml_dict(yaml_dict):

    def process_arg(arg):
        if isinstance(arg, list):
            return tuple(arg)
        return tuple([arg])

    # Get the names of the postprocesses to run
    if isinstance(yaml_dict["postprocess"], str):
        postprocess_names = yaml_dict["postprocess"].split(",")
    else:
        postprocess_names = yaml_dict["postprocess"]

    # Ensure that each name in postprocess_names exists in yaml_dict
    missing = [name for name in postprocess_names if name not in yaml_dict]
    if missing:
        raise ValueError(f"The following postprocess names are not defined in the YAML: {missing}")

    # Iterate through each postprocess
    out_postprocess = {}
    for name in postprocess_names:
        args = []; arg_values = []
        for arg in yaml_dict[name].get("combo_args", []):
            args.append(arg)
            arg_values.append(process_arg(yaml_dict[name]["combo_args"][arg]))

        arg_combos = list(it.product(*arg_values))

        for added_combo in yaml_dict[name].get("add_combo", []):
            combo = []
            for arg in args:
                if arg not in yaml_dict[name]["add_combo"][added_combo]:
                    raise ValueError(f"add_combo '{added_combo}' is missing key '{arg}'")
                combo.append(yaml_dict[name]["add_combo"][added_combo][arg])
            if tuple(combo) not in arg_combos:
                arg_combos.append(tuple(combo))

        for rm_combo in yaml_dict[name].get("ignore_combo", []):
            combo = []
            for arg in args:
                if arg not in yaml_dict[name]["ignore_combo"][rm_combo]:
                    raise ValueError(f"ignore_combo '{rm_combo}' is missing key '{arg}'")
                combo.append(yaml_dict[name]["ignore_combo"][rm_combo][arg]) 

            if tuple(combo) in arg_combos:
                arg_combos.remove(tuple(combo))

        postprocess_args = [
            {arg: val for arg, val in zip(args, combo)}
            for combo in arg_combos
        ]

        out_postprocess[name] = postprocess_args

    return out_postprocess

def create_df(postprocess_dict, existing_df=None):

    def get_dir(index):
        i = index + 1
        return f"{i:03d}"

    if existing_df is None:
        dfs = []
    else:
        dfs = [existing_df]

    arg_list = []
    for name, args in postprocess_dict.items():

        tmp = pd.DataFrame(args)
        arg_str = ",".join(tmp.columns)
        
        for arg in tmp.columns:
            if arg not in arg_list:
                arg_list.append(arg)

        tmp["name"] = name
        tmp["args"] = arg_str
        tmp["directory"] = None
        tmp["status"] = "new"

        dfs.append(tmp)

    for df in dfs:
        for arg in arg_list:
            if arg not in df.columns:
                df[arg] = None

    tmp = pd.concat(dfs)
    tmp = tmp[["name","directory","args","status"]+arg_list]

    dfs = []
    for name, name_df in tmp.groupby("name"):
        name_df = name_df.reset_index(drop=True).copy()
        name_df["directory"] = name_df.apply(lambda x: f"{name}/{get_dir(x.name)}" if x.directory is None else x.directory, axis=1)
        dfs.append(name_df)

    out_df = pd.concat(dfs)

    # Drop duplicates
    out_df = out_df.drop_duplicates(subset=["name"]+arg_list)
    
    return out_df

def postprocess_init(yaml_file):

    exp_args = load_yaml(yaml_file)

    exp_dir = Path(exp_args["experiment_directory"])

    tracking_file = exp_dir / "postprocess.tsv"

    if tracking_file.exists():
        tracking_df = pd.read_csv(tracking_file, sep="\t")
    else:
        tracking_df = None

    postprocess_runs = process_yaml_dict(exp_args)

    tracking_df = create_df(postprocess_runs, tracking_df)

    with open(tracking_file, "w") as outdf:
        tracking_df.to_csv(outdf, sep="\t", index=False)

    # Write out top-level post-processing
    postprocess_names = exp_args["postprocess"]
    for postprocess in postprocess_names:
        for arg in ["combo_args", "add_combo", "ignore_combo"]:
            if arg in exp_args[postprocess]:
                del exp_args[postprocess][arg]

    out_args = {i: exp_args[i] for i in postprocess_names}

    with open(exp_dir / 'postprocess.yaml', 'w') as outf:
        yaml.dump(out_args, outf)

def postprocess_commands(yaml_file):

    exp_args = load_yaml(yaml_file)

    exp_dir = Path(exp_args["experiment_directory"])

    with open(exp_dir / "yaml_files/yaml_files.txt", "r") as yfs:
        exp_list = [i.strip().replace(".yaml", "") for i in yfs.readlines()]
        exp_list = ["hehe", "haha"]

    tracking_file = exp_dir / "postprocess.tsv"

    if tracking_file.exists():
        tracking_df = pd.read_csv(tracking_file, sep="\t")
    else:
        # exit
        print("Need to init") # change to more informative msg

    to_do_df = tracking_df[tracking_df.status.isin(["new", "rerun"])]

    base_cmd = f"python blah blah"
    for _, row in to_do_df.iterrows():

        name = row["name"]
        
        cmd = [f"--set override_yaml={exp_dir}/postprocess.yaml",
               f"--set postprocess={name}"]

        for arg in row["args"].split(","):
            cmd.append(f"--set {name}.{arg}={row[arg]}")

        for run in exp_list:
            exp_cmd = [f"{run}/args.yaml"]
            print(" ".join([base_cmd] + exp_cmd + cmd))

    # TODO: write to a similar bash script as the actual experiment manager that includes all the commands

def postprocess_describe(yaml_file):
    pass



def main():
    import argparse

    parser = argparse.ArgumentParser(
        prog="postprocess_experiment.py",
        description="IBD-sims post-processing manager for experiments.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # init
    p_init = sub.add_parser("init", help="Initialize or re-initialize post-processing runs.")
    p_init.add_argument("yaml", help="Path to post-processing YAML")

    # commands
    p_cmd = sub.add_parser("commands", help="Print postprocess.py commands for simulations")
    p_cmd.add_argument("yaml", help="Path to post-processing YAML")

    # commands
    p_desc = sub.add_parser("describe", help="Describe the set of postprocessing runs.")
    p_desc.add_argument("yaml", help="Path to post-processing YAML")


    args = parser.parse_args()
    exp_dict = load_yaml(args.yaml)

    if args.command == "init":
        postprocess_init(args.yaml)

    if args.command == "commands":
        postprocess_commands(args.yaml)

if __name__ == "__main__":
    main()

    """
    Questions:
    
    1. How does the current code deal with --set hapne_ibd.filter=nan
    2. Can we add an option to use a yaml to override or add to an existing yaml? Like override_yaml is trying to do
    """