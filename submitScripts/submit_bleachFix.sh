#!/bin/bash
imagingDataDir=$1
behaviorDataDir=$2
matlabPath=$3
flyNum=$4

echo $imagingDataDir
echo $behaviorDataDir
echo $matlabPath

cwd=$(pwd)
parentdir="$(dirname "$cwd")"
traceFolder="$imagingDataDir"Yproj/""

newBleachBuffer=4

# # trim imaging data as if bleachBuffer in runs 2-end were nonzero
# $matlabPath -nodisplay -nodesktop -r "cd('../compilation/'); retroactivelyAddBleachBuffer $traceFolder $newBleachBuffer; exit"

# matlab realigns behavior traces to match imaging
$matlabPath -nodisplay -nodesktop -r "cd('../compilation/'); alignImagingAndBehaviorMultiImSingleBeh $parentdir $traceFolder $newBleachBuffer; exit"


# python smooths imaging, behavior, computes dFF and clustering
python ../behavior/smoothBehavior.py $traceFolder
python ../imaging/postProcessing/postProcessOne.py $traceFolder

