#!/bin/bash
#SBATCH -A axs
#SBATCH -J rs2
#SBATCH -t 2:00:00
#SBATCH -c 2
#SBATCH -o RS2.out -e RS2.err
#SBATCH --mem-per-cpu=4gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ess2129@columbia.edu

# sleep 3600

i=50
for f in /moto/axs/projects/rawdata/fly2/*.mat
do
	module load matlab/2018b
	# export MATLAB_PREFDIR= /moto/axs/users/code/jobs/
	echo "Launching a Matlab run"
	date
	let "i+=1"
	echo $i
	echo $f
	sbatch -J s$i	\
	-A axs	\
	-t 11:59:00	\
	-c 12	\
	--mem-per-cpu=8gb	\
	-o evanOut/s$i.out	\
	--wrap="matlab -nosplash -nodisplay -nodesktop -r 'registerMatfileDualColor_V5('\''$f'\'')' # > matoutfile"	
	sleep 30
done
# End of script
