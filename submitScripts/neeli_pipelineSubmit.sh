#!/bin/bash

# --- these are the only lines to edit ------------------------------
# -------------------------------------------------------------------
flyNum="fly3"
imagingDataDir="/Volumes/SCAPEdata1/finalData/2018_08_24_NsybNLS_odors/"$flyNum"/"
behaviorDataDir=""
matlabPath="/Applications/MATLAB_R2018b.app/bin/matlab"
# -------------------------------------------------------------------
# -------------------------------------------------------------------
bash submit_postProcessing_tmp.sh $imagingDataDir $behaviorDataDir $matlabPath $flyNum