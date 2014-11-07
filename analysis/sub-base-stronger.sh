#!/bin/sh
#PBS -l walltime=05:00:00
#PBS -M sean_anderson@sfu.ca
#PBS -l pmem=650mb
#PBS -l procs=1

module load R/3.1.1
cd $PBS_O_WORKDIR
R CMD BATCH --no-save 2-run-models-base-stronger-prior.R

