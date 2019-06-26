
function alignImagingAndBehaviorMultiTrial(codePath, traceFolder)
% In 'experiment' folder, expects to find concatenated imaging file F.mat.
% Also expects a folder called 'behavior' containing each of the mat files
% ouput from extractBehaviorManual.m, named flyX_runY.mat, and a similarly
% populated folder called info

addpath(genpath(codePath))

%traceFolder = [baseFolder,experiment];%[baseFolder,'rustyOut/',experiment];
imagingFile = [traceFolder,'F.mat']; %[traceFolder,'all.mat'];
load(imagingFile,'trialFlag');
behaviorDirectory = dir([traceFolder,'behavior/fly*']); %dir([baseFolder,'scapeBehavior/',experiment,'*.mat']);
extraSkippedFrames = 0; % DEPRECATED. needs to be vector if reintegrated
    
% get file order
tmp = zeros(size(behaviorDirectory));
for j=1:length(behaviorDirectory)
    tmp(j)=str2double(behaviorDirectory(j).name(end-4));
end
if sum(isnan(tmp))>0; error('Invalid File Order'); end
[runIds,fileOrder] = sort(tmp,'ascend');

alignedBehaviorTot = struct;

for i=1:length(behaviorDirectory)
    
    % point to behavior file (load below)
    behaviorFilename = [behaviorDirectory(fileOrder(i)).folder,'/',behaviorDirectory(fileOrder(i)).name];
    
    % load info file
    infoFile = [traceFolder,'info/',behaviorDirectory(fileOrder(i)).name(1:end-4),'_info.mat'];
    load(infoFile,'info');
    
    % grab length of correct piece of imaging data
    imRunLength = sum(trialFlag==runIds(i));
    
    % time
    bleachBuffer = 30*round(info.daq.scanRate);
    if info.daq.numberOfScans-bleachBuffer-extraSkippedFrames~=imRunLength; error('length of traces does not match expected'); end
    timeTot = (1:info.daq.numberOfScans)/info.daq.scanRate; %1/info.daq.scanRate:1/info.daq.scanRate:info.daq.scanLength;
    time = timeTot(bleachBuffer+extraSkippedFrames+1:end);
    %time = time(extraSkippedFrames+1:end);
    
    % behavior
    behaviorOpts.basicBehav = true; %false; % true => old defaults, ok if behav is const 30fps
    behaviorOpts.info = info;
    behaviorOpts.timeTot = timeTot;
    behaviorOpts.bleachBuffer = bleachBuffer;
    
    beh = formatBehaviorData(behaviorFilename, behaviorOpts);
    alignedBehaviorTot = concatenateBehaviorFiles(alignedBehaviorTot, beh, time);

end

bmat = matfile([traceFolder,'alignedBehavAndStim.mat'],'Writable',true);
bmat.alignedBehavior = alignedBehaviorTot;
bmat.time = alignedBehaviorTot.timeTot;



function behAligned = formatBehaviorData(behaviorFilename,behaviorOpts)

% takes as input a matfile made by extractBehaviorManual.m

beh = load(behaviorFilename);
beh.time.imOn = find(beh.traces.isImagingOn,1);
beh.time.imOff = find(beh.traces.isImagingOn,1,'last');
beh.smoothing = max(round( ((beh.time.imOff-beh.time.imOn)/behaviorOpts.info.daq.scanLength)/behaviorOpts.info.daq.scanRate ),1);    % smooth by behavior frames per imaging frame
beh.time.trueTime = (1:(beh.time.imOff-beh.time.imOn+1))/(beh.time.imOff-beh.time.imOn+1)*behaviorOpts.info.daq.scanLength;            % fixed?

beh.traces.legVarSmooth = smooth(beh.traces.legVar(beh.time.imOn:beh.time.imOff), beh.smoothing)';
beh.traces.stimSmooth = beh.traces.isStimOn(beh.time.imOn:beh.time.imOff);
beh.traces.drinkSmooth = beh.traces.isDrinking(beh.time.imOn:beh.time.imOff);

beh.traces.legVarSmoothAligned = interp1(beh.time.trueTime, beh.traces.legVarSmooth, behaviorOpts.timeTot);    % interpolate to align time
beh.traces.stimSmoothAligned = interp1(beh.time.trueTime, beh.traces.stimSmooth, behaviorOpts.timeTot);    % interpolate to align time
beh.traces.drinkSmoothAligned = interp1(beh.time.trueTime, beh.traces.drinkSmooth, behaviorOpts.timeTot);    % interpolate to align time

behAligned.legVar = beh.traces.legVarSmoothAligned(behaviorOpts.bleachBuffer+1:end);
behAligned.stim = beh.traces.stimSmoothAligned(behaviorOpts.bleachBuffer+1:end);
behAligned.drink = beh.traces.drinkSmoothAligned(behaviorOpts.bleachBuffer+1:end);



function alignedBehaviorTot = concatenateBehaviorFiles(alignedBehaviorTot, alignedBehavior, time)

% same as standalone function concatenateBehaviorTrials.m, but without
% loading and without loop

if ~isfield(alignedBehaviorTot,'legVar')
    alignedBehaviorTot.legVar = alignedBehavior.legVar; %behSmooth; %
    alignedBehaviorTot.stim = alignedBehavior.stim;
    alignedBehaviorTot.drink = alignedBehavior.drink;
    alignedBehaviorTot.timeTot = time;
else
    alignedBehaviorTot.legVar = [alignedBehaviorTot.legVar,alignedBehavior.legVar]; %behSmooth; %
    alignedBehaviorTot.stim = [alignedBehaviorTot.stim,alignedBehavior.stim];
    alignedBehaviorTot.drink = [alignedBehaviorTot.drink,alignedBehavior.drink];
    alignedBehaviorTot.timeTot = [alignedBehaviorTot.timeTot,alignedBehaviorTot.timeTot(end)+time];
end



