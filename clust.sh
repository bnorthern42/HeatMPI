#!/bin/bash

#SBATCH --job-name=bn4hsm
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:20:00
#SBATCH --mail-user=banorthern42@tntech.edu
#SBATCH --mail-type=ALL
#SBATCH --account=csc6740-001-2023f
#SBATCH --partition=batch-warp
#SBATCH --output=job_c8.txt
#SBATCH --mem=10G
#spack load /w22ygzf

spack unload --all
module load gnu9
module load openmpi4
cd /home/tntech.edu/banorthern42/heat2dmpi
make clean
make cart
mpirun -np 10 ./cartesian.out 2000 4000 10
