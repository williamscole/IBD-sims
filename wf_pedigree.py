import numpy as np
import sys
import argparse

def parse_arguments():
    # Create the parser
    parser = argparse.ArgumentParser(description='msprime simulation script')

    parser.add_argument("--constant_Ne", default=10_000, type=int)

    parser.add_argument("--samples", default=1000, type=int)

    parser.add_argument("--gens", default=25, type=int)

    parser.add_argument("--output", type=str)

    parser.add_argument("--monogamous", action="store_true")

    return parser.parse_args()

def add_di_parents(cur_gen, g, Ne, output, is_sample=False):
    n_samples = len(cur_gen)

    dads = np.random.randint(0, Ne // 2, n_samples)
    moms = np.random.randint(Ne // 2, Ne, n_samples)

    out_parents = set(dads) | set(moms)

    lines = [f"g{g-1}_{int(cur_gen[i])} g{g}_{int(dads[i])} g{g}_{int(moms[i])} {g-1} {1 if is_sample else 0} pop_0\n" 
             for i in range(n_samples)]

    output.writelines(lines)

    return list(out_parents)

def add_mono_parents(cur_gen, g, Ne, output, is_sample=False):
    n_samples = len(cur_gen)

    dads = np.random.randint(0, Ne // 2, n_samples)
    moms = (np.arange(Ne // 2) + (Ne//2))[dads]

    out_parents = set(dads) | set(moms)

    lines = [f"g{g-1}_{int(cur_gen[i])} g{g}_{int(dads[i])} g{g}_{int(moms[i])} {g-1} {1 if is_sample else 0} pop_0\n" 
             for i in range(n_samples)]

    output.writelines(lines)

    return list(out_parents)


def print_missing(cur_gen, g, output):
    for i in cur_gen:
        output.write(f"g{int(g)}_{int(i)} . . {int(g)} 0 pop_0\n")

def create_pedigree(samples: int,
                    demo_debug,
                    output: str,
                    gen: int = -1,
                    mating_func: int = "di",
                    seed: int = np.random.randint(1, 1000000),
                    population: str = "pop0"):

    cur_gen = np.arange(samples)

    output = open(f"{output}_WF.pedigree", "w")
    output.write("# id parent0 parent1 time is_sample population\n")

    for index, pop in enumerate(demo_debug.demography.populations):
        if pop.name == population:
            break

    for g in range(1, gen+1):

        np.random.seed(seed)
        seed += 10

        Ne = demo_debug.population_size_trajectory([g])[0,index]

        print(f"time {g}: {Ne} individuals")

        if mating_func == "di":
            cur_gen = add_di_parents(cur_gen, g, Ne, output, is_sample=g==1)

        elif mating_func == "mono":
            cur_gen = add_mono_parents(cur_gen, g, Ne, output, is_sample=g==1)

    print_missing(cur_gen, g, output)

    output.close()
        
    return np.arange(1, gen+1), demo_debug.population_size_trajectory(np.arange(1, gen+1))[:,index].astype(int)






if __name__ == "__main__":

    args = parse_arguments()

    if args.monogamous:
        monog(args)

    else:
        nonmonog(args)

