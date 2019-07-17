 #!/bin/bash

# --- these are the only lines to edit ------------------------------
# -------------------------------------------------------------------
exp="2019_06_28_Nsyb_NLS6s_walk/fly3"
dataDir="/Volumes/SCAPEdata1/finalData/"$exp"/"
saveDir="/Volumes/SCAPEdata1/figsAndMovies/"$exp"/"
matlabPath="/Applications/MATLAB_R2018b.app/bin/matlab"
# -------------------------------------------------------------------
# -------------------------------------------------------------------
bash submit_visualization_neeliTemp.sh $dataDir $saveDir $matlabPath $exp