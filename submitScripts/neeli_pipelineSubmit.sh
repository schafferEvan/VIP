#!/bin/bash

# --- these are the only lines to edit ------------------------------
# -------------------------------------------------------------------
flyNum="fly3"
imagingDataDir="/Volumes/SCAPEdata1/finalData/2019_06_28_Nsyb_NLS6s_walk/"$flyNum"/"
behaviorDataDir="/Volumes/SCAPEdata1/scapeBehavior/2019_06_28_Nsyb_NLS6s_walk/"
matlabPath="/Applications/MATLAB_R2018b.app/bin/matlab"
# -------------------------------------------------------------------
# -------------------------------------------------------------------
bash submit_postProcessing_tmp.sh $imagingDataDir $behaviorDataDir $matlabPath $flyNum