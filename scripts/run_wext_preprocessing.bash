#!/usr/bin/env bash

num_permutations=200
num_cores=$2
filename=$1
python ../wext/process_mutations.py \
    -m  $1 \
    -ct NA \
    -o  ${1}.data.json

python ../wext/compute_mutation_probabilities.py \
    -mf ${1}.data.json \
    -np $num_permutations \
    -nc $num_cores \
    -wf ${1}.weights.npy \
    -v  1
