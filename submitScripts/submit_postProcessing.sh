#!/bin/bash

cwd=$(pwd)
parentdir="$(dirname "$cwd")"
dataDir="/Volumes/SCAPEdata1/finalData/2019_05_01_57C10_NLS6f_Dilp5_GtACR/"

/Applications/MATLAB_R2018a.app/bin/matlab -nodisplay -nodesktop -r "cd('../imaging/postProcessing/'); compile_sumImage $parentdir $dataDir; exit"
