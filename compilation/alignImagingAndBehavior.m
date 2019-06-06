
function alignImagingAndBehavior
% -------------------------------------------------------------------------
% CHOOSE A DATASET
% -------------------------------------------------------------------------
addpath(genpath('~/Dropbox/_code/'))

experiment = '190422_f1/run4/'; %'2018_08_24_odor';
baseFolder = '/Users/evan/Dropbox/_sandbox/sourceExtraction/good/_runAndFeed/'; %'/Volumes/dataFast2/';%,'/fly3run1/'];

% -------------------------------------------------------------------------
traceFolder = [baseFolder,experiment];%[baseFolder,'rustyOut/',experiment];
infoFile = dir([traceFolder,'info/*info.mat']);
if length(infoFile)>1; infoFile=infoFile(end); end
load([infoFile.folder,'/',infoFile.name],'info');
imagingFile = [traceFolder,'F.mat']; %[traceFolder,'all.mat'];
scapeData = load(imagingFile);
behaviorDirectory = dir([traceFolder,'behavior.mat']); %dir([baseFolder,'scapeBehavior/',experiment,'*.mat']);
for i=1:length(behaviorDirectory); if behaviorDirectory(i).bytes>4096; break; end; end; behaviorDirectory = behaviorDirectory(i); % ignore .files


% check where to begin (if imaging traces contain multiple trials)
if exist([traceFolder,'info/skip/'],'dir')
    skipInfoFile = dir([traceFolder,'info/skip/*info.mat']);
    skipInfo = load([traceFolder,'info/skip/',skipInfoFile(2).name],'info');
    bleachBuffer = 30*round(skipInfo.info.daq.scanRate);
    framesToSkip = skipInfo.info.daq.numberOfScans-bleachBuffer;
else
    framesToSkip = 0;
end

% [dfF,F] = calculate_dff(scapeData.A, scapeData.C(:,framesToSkip+1:end), scapeData.YrA(:,framesToSkip+1:end), scapeData.b, scapeData.f(:,framesToSkip+1:end),0);
% %scapeData.Y = scapeData.C(:,framesToSkip+1:end) + scapeData.YrA(:,framesToSkip+1:end);
% scapeData.Y = F;
% scapeData.dfF = dfF;

% temp --
info.daq.numberOfScans = info.daq.numberOfScans;
%disp('subtracting 3729 frames for 0128 dataset')
extraSkippedFrames = 0;%3729; %=30 for 0226

% time
bleachBuffer = 30*round(info.daq.scanRate);
if info.daq.numberOfScans-bleachBuffer-extraSkippedFrames~=size(scapeData.F,2); error('length of traces does not match expected'); end
timeTot = (1:info.daq.numberOfScans)/info.daq.scanRate; %1/info.daq.scanRate:1/info.daq.scanRate:info.daq.scanLength;
time = timeTot(bleachBuffer+extraSkippedFrames+1:end);
%time = time(extraSkippedFrames+1:end);


% behavior
behaviorOpts.basicBehav = true; %false; % true => old defaults, ok if behav is const 30fps
behaviorOpts.info = info;
behaviorOpts.timeTot = timeTot;
behaviorOpts.bleachBuffer = bleachBuffer;

beh = formatBehaviorData(behaviorDirectory,behaviorOpts);
%alignedBehavior = beh.alignedBehavior(extraSkippedFrames+1:end);disp('subtracting 3729 frames for 0128 dataset')
%alignedStim = beh.alignedStim(extraSkippedFrames+1:end);disp('subtracting 3729 frames for 0128 dataset')

bmat = matfile([traceFolder,'alignedBehavAndStim.mat'],'Writable',true);
bmat.alignedBehavior = beh;
bmat.time = time;






function behAligned = formatBehaviorData(behaviorDirectory,behaviorOpts)

% takes as input a matfile made by extractBehaviorManual.m

beh = load([behaviorDirectory.folder,'/',behaviorDirectory.name]);
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


