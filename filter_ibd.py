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
                map_df = pd.read_csv(yargs["hapmap_chr1"].replace("chr1", f"chr{i}"))
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
        ordered_nodes = sorted([[added_nodes.index(i),i] for i in nodes], key=lambda x: x[0])
        to_remove = ordered_nodes[-1][1]
        nodes.discard(to_remove)
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
            degree_list = sorted([[sub.degree(i),i] for i in sub_nodes if i not in cur_nodes], key=lambda x: x[0])
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
    # print(f"Mean: {np.mean(ks)}, max: {max(ks)}, number: {len(out_nodes)}")

    return np.random.choice(out_nodes, target, replace=False) if len(out_nodes) > target else out_nodes


if __name__ == "__main__":

    ibd_df = pd.read_csv(sys.argv[1], delim_whitespace=True, header=None)
    n_samples = int(sys.argv[2])
    filtering = sys.argv[3] if len(sys.argv) > 3 else "none"

    if filtering == "none" or filtering == "null" or filtering == "":
        ibd_df.to_csv(sys.stdout, sep=" ", header=False, index=False)
        sys.exit()

    kinship_obj = Kinship(os.path.dirname(sys.argv[1]))

    rows = []
    for pair, pair_df in ibd_df.groupby([0, 2]):
        if pair_df[7].sum()>10:
            rows.append([*pair, pair_df[7].sum()])

    kinship = pd.DataFrame(rows, columns=["node1", "node2", "k"]).sort_values("k", ascending=False)
    kinship["related"] = kinship["k"].apply(kinship_obj.assign_related)
    related = kinship[kinship.related]

    all_nodes = {f"tsk_{i}" for i in range(n_samples)}

    if filtering == "related":

        nodes = subset_close(related, all_nodes, target=n_samples // 10)

    elif filtering == "unrelated":
        g = nx.Graph()

        g.add_nodes_from(list(all_nodes))

        g.add_edges_from(related[["node1", "node2"]].values)

        nodes = prune_relatives(g, target=n_samples // 10)

    ibd_df = ibd_df[ibd_df.apply(lambda x: x[0] in nodes and x[2] in nodes, axis=1)]

    ibd_df.to_csv(sys.stdout, sep=" ", header=False, index=False)



