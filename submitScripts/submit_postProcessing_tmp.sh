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

# matlab realigns behavior traces to match imaging
$matlabPath -nodisplay -nodesktop -r "cd('../compilation/'); alignImagingAndBehaviorMultiImSingleBeh $parentdir $traceFolder; exit"


# python smooths imaging, behavior, computes dFF and clustering
python ../behavior/smoothBehavior.py $traceFolder
python ../imaging/postProcessing/postProcessOne.py $traceFolder

