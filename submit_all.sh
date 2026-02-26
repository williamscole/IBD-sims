#!/bin/bash

for arg in "$@"; do
    bash pipeline.sh $arg
    echo "Submitted ${arg}!"
done
