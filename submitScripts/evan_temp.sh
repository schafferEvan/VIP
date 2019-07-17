#!/bin/bash

# --- these are the only lines to edit ------------------------------
# -------------------------------------------------------------------
flyNum="fly1"
# imagingDataDir="/Volumes/SCAPEdata1/finalData/2019_06_30_Nsyb_NLS6s_walk/"$flyNum"/"
imagingDataDir="/Volumes/dataFast/2019_06_30_Nsyb_NLS6s_walk/2019_06_30/"$flyNum"/"
behaviorDataDir="/Volumes/SCAPEdata1/scapeBehavior/2019_06_30_Nsyb_NLS6s_walk/"
matlabPath="/Applications/MATLAB_R2018a.app/bin/matlab"
# -------------------------------------------------------------------
# -------------------------------------------------------------------
bash temp.sh $imagingDataDir $behaviorDataDir $matlabPath $flyNum