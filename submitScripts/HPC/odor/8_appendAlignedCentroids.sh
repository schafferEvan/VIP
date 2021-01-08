#!/bin/bash
#SBATCH -A axs
#SBATCH -J centroids_8
#SBATCH -t 11:00:00
#SBATCH -c 12
#SBATCH -o nrun.out -e nrun.err
#SBATCH --mem-per-cpu=8gb

module load anaconda/3-5.3.1
python append_aligned_centroids.py
# End of script
