#!/bin/bash

# --- these are the only lines to edit ------------------------------
# -------------------------------------------------------------------
imagingDataDir="/Volumes/SCAPEdata1/finalData/2019_06_28_Nsyb_NLS6s_walk/fly3/"
behaviorDataDir="/Volumes/SCAPEdata1/scapeBehavior/2019_06_28_Nsyb_NLS6s_walk/"
matlabPath="/Applications/MATLAB_R2018b.app/bin/matlab"
# -------------------------------------------------------------------
# -------------------------------------------------------------------


cwd=$(pwd)
parentdir="$(dirname "$cwd")"

# matlab parses ROIs using watershed, then extracts timeseries
$matlabPath -nodisplay -nodesktop -r "cd('../imaging/postProcessing/'); compile_sumImage $parentdir $imagingDataDir; exit"


# copy behavior traces into subdirectory of imaging to aggregate final output
traceFolder="$imagingDataDir"/Yproj/""
cp -r info $traceFolder
behavTraceFolder="$traceFolder"/behavior/""
mkdir $behavTraceFolder
cp "$behaviorDataDir"*.mat"" $behavTraceFolder


# matlab realigns behavior traces to match imaging
$matlabPath -nodisplay -nodesktop -r "cd('../compilation/'); alignImagingAndBehaviorMultiImSingleBeh $parentdir $imagingDataDir; exit"


# python smooths imaging, behavior, computes dFF and clustering
python ../behavior/smoothBehavior.py $traceFolder
python ../imaging/postProcessing/postProcessAll.py $traceFolder

