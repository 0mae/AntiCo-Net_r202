#!/usr/bin/bash

# Number of CPUs must be specified
#export OMP_NUM_THREADS=1

cd /temp

Rscript --vanilla --slave /temp/scripts/r/count_og_displacement_2.R $1 $2 $3
