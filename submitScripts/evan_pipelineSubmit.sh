#!/bin/bash

# --- these are the only lines to edit ------------------------------
# -------------------------------------------------------------------
imagingDataDir="/Volumes/SCAPEdata1/finalData/2019_06_26_Nsyb_NLS6s_walk/fly1/"
behaviorDataDir="/Volumes/SCAPEdata1/scapeBehavior/2019_06_26_Nsyb_NLS6s_walk/"
matlabPath="/Applications/MATLAB_R2018a.app/bin/matlab"
# -------------------------------------------------------------------
# -------------------------------------------------------------------
bash submit_postProcessing.sh $imagingDataDir $behaviorDataDir $matlabPath