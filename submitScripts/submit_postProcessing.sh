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

#make behavior traces
ballthresh=0.7 		# pixel threshold for ball vs not ball (quantile of blurred image)
indicatorHeight=10  # number of rows at top of image from which to measure indicator signal
indicatorStart=1
$matlabPath -nodisplay -nodesktop -r "cd('../behavior/'); extractBehaviorAuto $behaviorDataDir $ballthresh $indicatorHeight $indicatorStart; exit"


# copy behavior traces into subdirectory of imaging to aggregate final output
cp -r "$imagingDataDir"info"" $traceFolder
behavTraceFolder="$traceFolder"behavior/""
mkdir $behavTraceFolder
cp "$behaviorDataDir$flyNum"*.mat"" $behavTraceFolder
dlcFolder="$behaviorDataDir"dlc/""
cp "$dlcFolder$flyNum"*.csv"" $behavTraceFolder



# matlab realigns behavior traces to match imaging (use this version for single behav file and version below for multiple)
$matlabPath -nodisplay -nodesktop -r "cd('../compilation/'); alignImagingAndBehaviorMultiImSingleBeh $parentdir $traceFolder; exit"

# # for experiments with multiple behavior files, align as follows instead of above
# $matlabPath -nodisplay -nodesktop -r "cd('../compilation/'); alignImagingAndBehaviorMultiTrial $parentdir $traceFolder; exit"


# # IF NECESSARY, throw away trials from processed files (i.e. if early data from a fly is good but late isn't)
# trialsToSkip = '4:8' # use matlab array syntax but pass as string
# $matlabPath -nodisplay -nodesktop -r "cd('../compilation/'); retroactivelyIgnoreTrials $traceFolder $trialsToSkip; exit"


# generate aggregated RAW mat files to share with collaborators
fromGreen='False'
python3 ../compilation/compileRaw.py $traceFolder $expID $traceFolder $fromGreen

# # python smooths imaging, behavior, computes dFF and clustering
python3 ../imaging/postProcessing/postProcessOne.py $traceFolder $expID 

