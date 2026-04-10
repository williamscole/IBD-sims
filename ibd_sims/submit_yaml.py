#!/usr/bin/env python3

import subprocess
import time
import argparse
import os
from typing import Set
import sys
import yaml
import datetime
import shutil



def write_job_list(path):
    args = yaml.safe_load(open(f"{path}/args.yaml"))

    print(args)

    out = []

    jobid = 1
    for i in range(args["iter"]):
        for c in range(args["end_chr"]-1):
            out.append(f"{jobid} {i+1} {c+2}")
            jobid += 1

    with open(f"{path}/array_list.txt", "w") as al:
        al.write("\n".join(out))


if __name__ == "__main__":
    write_job_list(sys.argv[1])