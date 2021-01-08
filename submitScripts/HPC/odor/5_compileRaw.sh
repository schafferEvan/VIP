#!/bin/bash
#SBATCH -A axs
#SBATCH -J compileRaw_5
#SBATCH -t 11:00:00
#SBATCH -c 12
#SBATCH -o nrun.out -e nrun.err
#SBATCH --mem-per-cpu=8gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nm2786@columbia.edu

module load anaconda/3-5.3.1
python compileRaw.py

# End of script
