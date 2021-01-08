#!/bin/bash
#SBATCH -A axs
#SBATCH -J align_4
#SBATCH -t 11:00:00
#SBATCH -c 12
#SBATCH -o nrun.out -e nrun.err
#SBATCH --mem-per-cpu=8gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nm2786@columbia.edu

cp -r "/moto/axs/projects/registered/2018_08_24_fly2/fly2/info/" "/moto/axs/projects/registered/2018_08_24_fly2/fly2/Yproj/"
behavTraceFolder="/moto/axs/projects/registered/2018_08_24_fly2/fly2/Yproj/behavior/"
mkdir $behavTraceFolder
cp "/moto/axs/projects/behavior/2018_08_24_fly2/fly2/"*.mat"" $behavTraceFolder

module load matlab/2018b
matlab -nosplash -nodisplay -nodesktop -r alignImagingAndBehaviorMultiTrial

# End of script
