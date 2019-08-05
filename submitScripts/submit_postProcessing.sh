#!/bin/bash
imagingDataDir=$1
behaviorDataDir=$2
matlabPath=$3
flyNum=$4
expID=$5

echo $imagingDataDir
echo $behaviorDataDir
echo $matlabPath

cwd=$(pwd)
parentdir="$(dirname "$cwd")"
traceFolder="$imagingDataDir"Yproj/""


# matlab parses ROIs using watershed
$matlabPath -nodisplay -nodesktop -r "cd('../imaging/postProcessing/'); compile_sumImage $parentdir $imagingDataDir; exit"

# matlab extracts time series
$matlabPath -nodisplay -nodesktop -r "cd('../imaging/postProcessing/'); extract_F_from_conComp $parentdir $imagingDataDir; exit"

# make behavior traces
#$matlabPath -nodisplay -nodesktop -r "cd('../behavior/'); extractBehaviorAuto $behaviorDataDir; exit"


# copy behavior traces into subdirectory of imaging to aggregate final output
cp -r "$imagingDataDir"info"" $traceFolder
behavTraceFolder="$traceFolder"behavior/""
mkdir $behavTraceFolder
cp "$behaviorDataDir$flyNum"*.mat"" $behavTraceFolder



# # trim imaging data as if bleachBuffer in runs 2-end were nonzero
newBleachBuffer=7
$matlabPath -nodisplay -nodesktop -r "cd('../compilation/'); retroactivelyAddBleachBuffer $traceFolder $newBleachBuffer; exit"

# matlab realigns behavior traces to match imaging
$matlabPath -nodisplay -nodesktop -r "cd('../compilation/'); alignImagingAndBehaviorMultiImSingleBeh $parentdir $traceFolder $newBleachBuffer; exit"

# generate aggregated RAW mat files to share with collaborators
fromGreen='False'
python ../compilation/compileRaw.py $traceFolder $expID $traceFolder $fromGreen


# python smooths imaging, behavior, computes dFF and clustering
python ../behavior/smoothBehavior.py $traceFolder
python ../imaging/postProcessing/postProcessOne.py $traceFolder


# generate aggregated npz files to share with collaborators
fromGreen='False'
python ../compilation/compileFinalSummary.py $traceFolder $expID $traceFolder $fromGreen