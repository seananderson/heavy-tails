#!/bin/sh
#PBS -l walltime=01:00:00
#PBS -M sean_anderson@sfu.ca
#PBS -l pmem=2200mb
#PBS -l procs=1

module load R/3.1.1
cd $PBS_O_WORKDIR
R CMD BATCH --no-save 5.8-stan-beta-modelling.R
