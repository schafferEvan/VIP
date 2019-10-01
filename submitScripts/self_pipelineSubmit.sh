 #!/bin/bash

# --- these are the only lines to edit ------------------------------
# -------------------------------------------------------------------
flyNum="fly3"
expName="2019_08_13_Nsyb_NLS6s_walk/"
imagingDataDir="/Volumes/SCAPEdata1/finalData/"$expName$flyNum"/"
behaviorDataDir="/Volumes/SCAPEdata1/scapeBehavior/"$expName""
matlabPath="/Applications/MATLAB_R2018b.app/bin/matlab"
datasharePath="/Volumes/dropboxdrive/Dropbox/00_AcademicFiles/12_ColumbiaDocumentsG6/00_AxelLab/01_Data/flygenvectors/datashare/SynologyDrive/"
# -------------------------------------------------------------------
# -------------------------------------------------------------------
bash submit_postProcessing.sh $imagingDataDir $behaviorDataDir $matlabPath $flyNum $expName $datasharePath