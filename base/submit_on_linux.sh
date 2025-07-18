#!/bin/bash -l
#SBATCH -o testMatCell.out
#SBATCH --qos=normal
#SBATCH -J TestMatCell
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH -t 24:00:00


module purge
module add matlab/r2023a


matlab -nodesktop -nodisplay -nosplash -r "run('../startup.m'); testMatCell; exit"
