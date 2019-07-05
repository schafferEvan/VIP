#!/bin/bash
imagingDataDir=$1
behaviorDataDir=$2
matlabPath=$3

echo $imagingDataDir
echo $behaviorDataDir
echo $matlabPath

cwd=$(pwd)
parentdir="$(dirname "$cwd")"

# matlab parses ROIs using watershed
$matlabPath -nodisplay -nodesktop -r "cd('../imaging/postProcessing/'); compile_sumImage $parentdir $imagingDataDir; exit"

# matlab extracts time series
$matlabPath -nodisplay -nodesktop -r "cd('../imaging/postProcessing/'); extract_F_from_conComp $parentdir $imagingDataDir; exit"

# make behavior traces
$matlabPath -nodisplay -nodesktop -r "cd('../imaging/postProcessing/'); extractBehaviorAuto $behaviorDataDir; exit"

# copy behavior traces into subdirectory of imaging to aggregate final output
traceFolder="$imagingDataDir"Yproj/""
cp -r info $traceFolder
behavTraceFolder="$traceFolder"/behavior/""
mkdir $behavTraceFolder
cp "$behaviorDataDir"*.mat"" $behavTraceFolder


# matlab realigns behavior traces to match imaging
$matlabPath -nodisplay -nodesktop -r "cd('../compilation/'); alignImagingAndBehaviorMultiImSingleBeh $parentdir $imagingDataDir; exit"


# python smooths imaging, behavior, computes dFF and clustering
python ../behavior/smoothBehavior.py $traceFolder
python ../imaging/postProcessing/postProcessAll.py $traceFolder

