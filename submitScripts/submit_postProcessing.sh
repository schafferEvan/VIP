#!/bin/bash

# --- these are the only lines to edit ------------------------------
# -------------------------------------------------------------------
dataDir="/Volumes/SCAPEdata1/finalData/2019_06_21_Nsyb_NLS6s_walk/fly2/"
matlabPath="/Applications/MATLAB_R2018a.app/bin/matlab"
# -------------------------------------------------------------------
# -------------------------------------------------------------------

cwd=$(pwd)
parentdir="$(dirname "$cwd")"
$matlabPath -nodisplay -nodesktop -r "cd('../imaging/postProcessing/'); compile_sumImage $parentdir $dataDir; exit"
python ../imaging/postProcessing/postProcessAll.py $dataDir