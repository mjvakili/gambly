#!/bin/bash

#PBS -l nodes=2:ppn=20,mem=400GB
#PBS -l walltime=40:00:00
#PBS -N Mr20-hod-group
#PBS -M mv1003@nyu.edu,mohammadjavad.vakili@gmail.com
#PBS -m abe
#PBS -j oe

module purge
module load mpi4py/openmpi/intel/1.3.1
module load numpy/intel/1.10.1
module load corrfunc/intel/1.0.0

cd /scratch/mv1003/projects/gambly/run

mpiexec -np 40 python /scratch/mv1003/projects/group20/gambly/code/adhoc_mcmc_group.py 140 30000 20.

# leave a blank line at the end

