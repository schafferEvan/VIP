 #!/bin/bash

# --- these are the only lines to edit ------------------------------
# -------------------------------------------------------------------
flyNum="fly2"
imagingDataDir="/Volumes/SCAPEdata1/finalData/2019_07_01_Nsyb_NLS6s_walk/"$flyNum"/"
behaviorDataDir="/Volumes/SCAPEdata1/scapeBehavior/2019_07_01_Nsyb_NLS6s_walk/"
matlabPath="/Applications/MATLAB_R2018a.app/bin/matlab"
# -------------------------------------------------------------------
# -------------------------------------------------------------------
bash submit_postProcessing.sh $imagingDataDir $behaviorDataDir $matlabPath $flyNum