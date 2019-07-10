 #!/bin/bash

# --- these are the only lines to edit ------------------------------
# -------------------------------------------------------------------
exp="2019_07_01_Nsyb_NLS6s_walk/fly1/"
dataDir="/Volumes/SCAPEdata1/finalData/"$exp""
saveDir="/Volumes/SCAPEdata1/figsAndMovies/"$exp""
matlabPath="/Applications/MATLAB_R2018a.app/bin/matlab"
# -------------------------------------------------------------------
# -------------------------------------------------------------------
bash submit_visualization.sh $dataDir $saveDir $matlabPath