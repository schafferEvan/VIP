#!/bin/bash
# NOTE: this can be run as a bash script, but it is primarily meant as a directory 
# of commands to dunk into a personalized submit script on an as-need basis

dataDir=$1
saveDir=$2
matlabPath=$3
expID=$4


echo $dataDir
echo $saveDir
echo $matlabPath

cwd=$(pwd)
parentdir="$(dirname "$cwd")"

# make quick temporally downsampled movie
$matlabPath -nodisplay -nodesktop -r "cd('../visualization/'); make_quickMovie $parentdir $dataDir $saveDir; exit"


# # make plots of traces of active cells with spatial footprints and behavior traces
# #fromGreenCC = false; 	 % use ROIs parsed from green image (false -> use red)
# #showBallVar = true;     % show motion energy of ball extracted as pix var
# #showDrink = false;      % show trace of bubble and drinking
# #showDLC = true;         % show output of deeplabcut
# fromGreenCC=0
# showBallVar=1
# showDrink=0
# showDLC=1
# $matlabPath -nodisplay -nodesktop -r "cd('../visualization/'); plotActiveNeurons $parentdir $dataDir $expID $saveDir $fromGreenCC $showBallVar $showDrink $showDLC; exit"


# make full movie (all frames)
# $matlabPath -nodisplay -nodesktop -r "cd('../visualization/'); make_rawScaled_movies $parentdir $dataDir $saveDir; exit"

# make component movie
# (work in progress)