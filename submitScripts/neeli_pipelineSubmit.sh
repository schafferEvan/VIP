 #!/bin/bash

# --- these are the only lines to edit ------------------------------
# -------------------------------------------------------------------
flyNum="fly3"
expDate="2019_08_20_Nsyb_NLS6s_walk/"
expID="$expDate$flyNum"
imagingDataDir="/Volumes/SCAPEdata1/finalData/"$expID"/"
behaviorDataDir="/Volumes/SCAPEdata1/scapeBehavior/"$expDate""
matlabPath="/Applications/MATLAB_R2018b.app/bin/matlab"
# -------------------------------------------------------------------
# -------------------------------------------------------------------
bash submit_postProcessing.sh $imagingDataDir $behaviorDataDir $matlabPath $flyNum $expID