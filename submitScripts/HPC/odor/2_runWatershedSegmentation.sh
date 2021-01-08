#!/bin/bash
#SBATCH -A axs
#SBATCH -J segment_2
#SBATCH -t 11:00:00
#SBATCH -c 12
#SBATCH -o nrun.out -e nrun.err
#SBATCH --mem-per-cpu=8gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nm2786@columbia.edu

module load matlab/2018b
matlab -nosplash -nodisplay -nodesktop -r run_watershed_segmentation

# End of script
