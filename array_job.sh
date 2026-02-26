#!/bin/bash

path=${1}
index=${2}
iter=${3}
end_chr=${4}

iter=$(( ($index-1)/${end_chr} + 1 ))
chrom=$(( ($index-1)%${end_chr} + 1 ))

echo "Running iteration $iter, chromosome $chrom"

python simulations.py ${path} ${iter} ${chrom}