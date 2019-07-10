#!/bin/bash
# NOTE: this can be run as a bash script, but it is primarily meant as a directory 
# of commands to dunk into a personalized submit script on an as-need basis

dataDir=$1
saveDir=$2
matlabPath=$3

echo $dataDir
echo $saveDir
echo $matlabPath

cwd=$(pwd)
parentdir="$(dirname "$cwd")"

# make quick temporally downsampled movie
$matlabPath -nodisplay -nodesktop -r "cd('../visualization/'); make_quickMovie $parentdir $dataDir $saveDir; exit"
