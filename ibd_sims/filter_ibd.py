import pandas as pd
import sys
import networkx as nx
import pickle as pkl
import itertools as it
import numpy as np
import yaml
import os


class Kinship:
    def __init__(self, path):
        yargs = yaml.safe_load(open(f"{path}/args.yaml"))

        end_chr = yargs["end_chr"]

        if end_chr == 30:
            self.tot = 3000
        elif end_chr == 22:
            self.tot = 0
            for i in range(1, 23):
                map_df = pd.read_csv(yargs["hapmap_chr1"].replace("chr1", f"chr{i}"), sep="\\s+")
                self.tot += map_df.iloc[-1]["Map(cM)"] - map_df.iloc[0]["Map(cM)"]
        elif end_chr == 2:
            self.tot = 200
        elif end_chr == 1:
            self.tot = 100

        first = (2**-.5*self.tot, self.tot)
        second = (2**-1.5*self.tot, 2**-.5*self.tot)
        third = (2**-2.5*self.tot, 2**-1.5*self.tot)

        self.degrees = {
            "1st": first,
            "2nd": second,
            "3rd": third
        }

    def assign_related(self, tot_ibd, related=["1st", "2nd", "3rd"]):
        degree = "UN"
        for d, (lw, up) in self.degrees.items():
            if lw < tot_ibd <= up:
                degree = d
        return degree in related


def subset_close(related_df, node_set, target=1000):
    nodes = set()
    added_nodes = list()
    for row in related_df.itertuples():
        if len(nodes) >= target:
            break
        nodes.add(row.node1)
        nodes.add(row.node2)
        added_nodes.extend([row.node1, row.node2])
    if len(nodes) == target + 1:
        nodes_sorted = sorted([[added_nodes.index(i), i] for i in nodes], key=lambda x: x[0])
        nodes = {n[1] for n in nodes_sorted[:-1]}
    elif len(nodes) < target:
        n_to_add = target - len(nodes)
        to_add = node_set - nodes
        nodes |= set(np.random.choice(list(to_add), n_to_add, replace=False))

    return nodes


def prune_relatives(g, target=1000, out_nodes=[]):
    cur_l = len(out_nodes)

    for sub_nodes in nx.connected_components(g):
        if len(sub_nodes) == 1:
            out_nodes.extend(sub_nodes)
            cur_l += 1
            continue

        sub = g.subgraph(sub_nodes)
        cur_nodes = []

        while True:
            degree_list = sorted([[sub.degree(i), i] for i in sub_nodes if i not in cur_nodes], key=lambda x: x[0])
            add = False

            for d, node1 in degree_list:
                add = True
                for node2 in cur_nodes:
                    if sub.has_edge(node1, node2):
                        add = False
                        break
                if add:
                    cur_nodes.append(node1)
                    cur_l += 1
                    break

            if add == False or cur_l >= target:
                break

        out_nodes.extend(cur_nodes)

    ks = []
    for i, j in it.combinations(out_nodes, r=2):
        k = g.get_edge_data(i, j, {"weight": 0})["weight"]
        ks.append(k)

    return np.random.choice(out_nodes, target, replace=False) if len(out_nodes) > target else out_nodes


# --- Node file helpers ---

def _node_file_path(ibd_path, label):
    """Return path to the node list file for a given label, next to the IBD file.

    e.g. /path/to/iter1.ibd.gz -> /path/to/iter1_unrelated.txt
    """
    base = ibd_path.replace(".ibd.gz", "").replace(".ibd", "")
    return f"{base}_{label}.txt"


def _write_node_file(nodes, ibd_path, label):
    path = _node_file_path(ibd_path, label)
    with open(path, "w") as f:
        for node in sorted(nodes):
            f.write(f"{node}\n")
    print(f"Wrote {len(nodes)} nodes to {path}")


def _read_node_file(ibd_path, label):
    """Load a node file if it exists, returning a set of node IDs or None."""
    path = _node_file_path(ibd_path, label)
    if os.path.exists(path):
        with open(path) as f:
            return set(line.strip() for line in f if line.strip())
    return None


# --- Public helpers for node resolution ---

def _normalize_filtering(filtering):
    """Normalize a filtering value to one of: 'none', 'random', 'related', 'unrelated'.

    'none' (or null/None/"") means no filtering — all samples are kept.
    'random' means randomly subsample to n_samples // 10.
    """
    if filtering == "null" or filtering == "" or filtering is None:
        return "none"
    if filtering not in ("none", "random", "related", "unrelated"):
        raise ValueError(f"Unknown filtering mode: {filtering!r}")
    return filtering


def get_node_file_path(ibd_path, filtering):
    """Return the path to the cached node file for a given filtering mode.

    e.g. /path/to/iter1.ibd.gz + "unrelated" -> /path/to/iter1_unrelated.txt
    """
    filtering = _normalize_filtering(filtering)
    return _node_file_path(ibd_path, filtering)


def get_nodes(ibd_path, n_samples, filtering):
    """Return the set of node IDs for a given filtering mode.

    Modes:
      'none'       — return all n_samples nodes (no filtering at all)
      'random'     — randomly subsample to n_samples // 10
      'related'    — subset enriched for close relatives (n_samples // 10)
      'unrelated'  — subset pruned of close relatives (n_samples // 10)

    For 'random', 'related', 'unrelated': uses cached node files if available,
    otherwise computes and caches them.
    """
    filtering = _normalize_filtering(filtering)

    # No filtering — return everybody
    if filtering == "none":
        return {f"tsk_{i}" for i in range(n_samples)}

    # For the subsetting modes, check cache first
    nodes = _read_node_file(ibd_path, filtering)

    if nodes is not None:
        print(f"Loaded {len(nodes)} {filtering} nodes from cache.")
        return nodes

    print(f"No cached node file found for '{filtering}', computing fresh.")
    all_nodes = {f"tsk_{i}" for i in range(n_samples)}
    target = n_samples // 10

    if filtering == "random":
        nodes = set(np.random.choice(list(all_nodes), target, replace=False))

    else:
        ibd_df_full = pd.read_csv(ibd_path, sep="\\s+", header=None)
        kinship_obj = Kinship(os.path.dirname(ibd_path))

        rows = []
        for pair, pair_df in ibd_df_full.groupby([0, 2]):
            if pair_df[7].sum() > 10:
                rows.append([*pair, pair_df[7].sum()])

        kinship = pd.DataFrame(rows, columns=["node1", "node2", "k"]).sort_values("k", ascending=False)
        kinship["related"] = kinship["k"].apply(kinship_obj.assign_related)
        related = kinship[kinship.related]

        if filtering == "related":
            nodes = subset_close(related, all_nodes, target=target)

        elif filtering == "unrelated":
            g = nx.Graph()
            g.add_nodes_from(list(all_nodes))
            g.add_edges_from(related[["node1", "node2"]].values)
            nodes = set(prune_relatives(g, target=target))

    _write_node_file(nodes, ibd_path, filtering)
    return nodes


# --- Main public functions ---

def write_samples(ibd_path, n_samples):
    """Compute and write node files for all three filtering modes.

    Creates iter{i}_random.txt, iter{i}_related.txt, iter{i}_unrelated.txt
    next to the IBD file. Safe to call multiple times — will overwrite existing files.
    """
    ibd_df = pd.read_csv(ibd_path, sep="\\s+", header=None)
    all_nodes = {f"tsk_{i}" for i in range(n_samples)}
    target = n_samples // 10

    # random — no kinship needed
    random_nodes = set(np.random.choice(list(all_nodes), target, replace=False))
    _write_node_file(random_nodes, ibd_path, "random")

    # related + unrelated both need kinship
    kinship_obj = Kinship(os.path.dirname(ibd_path))

    rows = []
    for pair, pair_df in ibd_df.groupby([0, 2]):
        if pair_df[7].sum() > 10:
            rows.append([*pair, pair_df[7].sum()])

    kinship = pd.DataFrame(rows, columns=["node1", "node2", "k"]).sort_values("k", ascending=False)
    kinship["related"] = kinship["k"].apply(kinship_obj.assign_related)
    related = kinship[kinship.related]

    related_nodes = subset_close(related, all_nodes, target=target)
    _write_node_file(related_nodes, ibd_path, "related")

    g = nx.Graph()
    g.add_nodes_from(list(all_nodes))
    g.add_edges_from(related[["node1", "node2"]].values)
    unrelated_nodes = set(prune_relatives(g, target=target))
    _write_node_file(unrelated_nodes, ibd_path, "unrelated")


def filter_ibd(ibd_path, n_samples, output, filtering="none"):

    filtering = _normalize_filtering(filtering)

    if filtering == "none":
        # No filtering — copy the IBD file as-is
        ibd_df = pd.read_csv(ibd_path, sep="\\s+", header=None)
        ibd_df.to_csv(output, sep=" ", header=False, index=False)
        return

    nodes = get_nodes(ibd_path, n_samples, filtering)

    ibd_df = pd.read_csv(ibd_path, sep="\\s+", header=None)
    ibd_df = ibd_df[ibd_df.apply(lambda x: x[0] in nodes and x[2] in nodes, axis=1)]
    ibd_df.to_csv(output, sep=" ", header=False, index=False)


if __name__ == "__main__":

    ibd_path = sys.argv[1]
    n_samples = int(sys.argv[2])
    filtering = sys.argv[3] if len(sys.argv) > 3 else "none"

    filter_ibd(ibd_path, n_samples, sys.stdout, filtering)