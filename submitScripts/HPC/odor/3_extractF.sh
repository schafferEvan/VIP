#!/bin/bash
#SBATCH -A axs
#SBATCH -J extractF_3
#SBATCH -t 11:00:00
#SBATCH -c 12
#SBATCH -o nrun.out -e nrun.err
#SBATCH --mem-per-cpu=8gb

module load matlab/2018b
matlab -nosplash -nodisplay -nodesktop -r extract_F_from_conComp

# End of script
