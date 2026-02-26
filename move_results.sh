#!/bin/bash

newdir=${1}

for arg in "${@:2}"; do
    mkdir ${arg}/${newdir}
    mv ${arg}/*.ne ${arg}/${newdir}
    mv ${arg}/*.pedigree ${arg}/${newdir}
    mv ${arg}/*.npy ${arg}/${newdir}
    mv ${arg}/*.ibd.gz ${arg}/${newdir}
    mv ${arg}/*.tmrca.gz ${arg}/${newdir}
    mv ${arg}/*.map ${arg}/${newdir}
    echo "Submitted ${arg}!"
done