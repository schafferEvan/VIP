
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
repeatedBleachBuffer = 0; % this is now handled in postprocessing

    
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
        
    % ignore first 30s of beginning of experiment, but not subsequent runs
    if i==1; bleachBuffer = 30*round(info.daq.scanRate);
    else;    bleachBuffer = repeatedBleachBuffer*round(info.daq.scanRate);
    end
    
    % check that extracted imaging data has correct number of frames
    disp({'# scans: ',num2str(info.daq.numberOfScans);...
        'bleachBuffer: ',num2str(bleachBuffer);...
        'imRunLength: ',num2str(imRunLength)})
    if info.daq.numberOfScans-bleachBuffer~=imRunLength
        s=strfind(traceFolder,'/');
        exp_folder = traceFolder(1:s(end-1));
        [~, trialOrder, ~, ~, frameNum] = sortExperimentDirectory(exp_folder,'reg');
        bleachBuffer = (frameNum(trialOrder(1)))-1; %str2double
        disp({'# scans: ',num2str(info.daq.numberOfScans);...
        'bleachBuffer_adjusted: ',num2str(bleachBuffer);...
        'imRunLength: ',num2str(imRunLength)})
        if (info.daq.numberOfScans-bleachBuffer==imRunLength)
            warning('bleach buffer inferred from matfile name');
        else
            error('length of traces does not match expected'); 
        end
    end
    
    % populate imaging time vector
    timeTot = (1:info.daq.numberOfScans)/info.daq.scanRate;
    timeTmp = timeTot(bleachBuffer+1:end);

    %     if info.daq.numberOfScans-bleachBuffer-extraSkippedFrames~=imRunLength; error('length of traces does not match expected'); end
    %     timeTot = (1:info.daq.numberOfScans)/info.daq.scanRate; %1/info.daq.scanRate:1/info.daq.scanRate:info.daq.scanLength;
    %     time = timeTot(bleachBuffer+extraSkippedFrames+1:end);
    %     %time = time(extraSkippedFrames+1:end);
    
    % behavior
    behRaw = load(behaviorFilename);
    behaviorOpts.parseStruct = getBehParsing(double(behRaw.traces.isImagingOn), traceFolder);
    
    % check for deeplabcut output and load if available
    dlcname = dir([traceFolder,'behavior/*DeepCut*.csv']);
    if ~isempty(dlcname)
        behRaw.nDLCpts = 8; disp('update dlcRead to allow for variable number of points')
        behRaw.dlcData = dlcRead([traceFolder, 'behavior/', dlcname.name],behRaw.nDLCpts);
    else
        behRaw.nDLCpts = 0;
        behRaw.dlcData = [];
    end
    
    behaviorOpts.basicBehav = true; %false; % true => old defaults, ok if behav is const 30fps
    behaviorOpts.info = info;
    behaviorOpts.timeTot = timeTot;
    behaviorOpts.bleachBuffer = bleachBuffer;
    
    beh = formatBehaviorData(behRaw, behaviorOpts, 1); % iter=1, not i
    disp('Behavior file added')
    alignedBehaviorTot = concatenateBehaviorFiles(alignedBehaviorTot, beh, timeTmp);
    %     beh = formatBehaviorData(behaviorFilename, behaviorOpts);
    %     alignedBehaviorTot = concatenateBehaviorFiles(alignedBehaviorTot, beh, time);

end

bmat = matfile([traceFolder,'alignedBehavAndStim.mat'],'Writable',true);
bmat.alignedBehavior = alignedBehaviorTot;
bmat.time = alignedBehaviorTot.timeTot;
disp('--- COMPLETED ALIGNMENT OF BEHAVIOR AND IMAGING ---')



function behAligned = formatBehaviorData(beh, behaviorOpts, iter)

% takes as input a matfile made by extractBehaviorManual.m or extractBehaviorAuto.m

beh.time.imOn = behaviorOpts.parseStruct.starts(iter);
beh.time.imOff = behaviorOpts.parseStruct.stops(iter);
beh.smoothing = max(round( ((beh.time.imOff-beh.time.imOn)/behaviorOpts.info.daq.scanLength)/behaviorOpts.info.daq.scanRate ),1);    % smooth by behavior frames per imaging frame
beh.time.trueTime = (1:(beh.time.imOff-beh.time.imOn+1))/(beh.time.imOff-beh.time.imOn+1)*behaviorOpts.info.daq.scanLength;            % fixed?

% smooth, align, and trim traces from simple behavior extraction
beh.traces.legVarSmooth = smooth(beh.traces.legVar(beh.time.imOn:beh.time.imOff), beh.smoothing)';
beh.traces.stimSmooth = beh.traces.isStimOn(beh.time.imOn:beh.time.imOff);
beh.traces.drinkSmooth = beh.traces.isDrinking(beh.time.imOn:beh.time.imOff);

beh.traces.legVarSmoothAligned = interp1(beh.time.trueTime, beh.traces.legVarSmooth, behaviorOpts.timeTot,'nearest','extrap');    % interpolate to align time
beh.traces.stimSmoothAligned = interp1(beh.time.trueTime, beh.traces.stimSmooth, behaviorOpts.timeTot,'nearest','extrap');    % interpolate to align time
beh.traces.drinkSmoothAligned = interp1(beh.time.trueTime, beh.traces.drinkSmooth, behaviorOpts.timeTot,'nearest','extrap');    % interpolate to align time

behAligned.legVar = beh.traces.legVarSmoothAligned(behaviorOpts.bleachBuffer+1:end);
behAligned.stim = beh.traces.stimSmoothAligned(behaviorOpts.bleachBuffer+1:end);
behAligned.drink = beh.traces.drinkSmoothAligned(behaviorOpts.bleachBuffer+1:end);


% smooth, align, and trim traces from deeplabcut behavior extraction
beh.dlcAligned.tmp.smooth = zeros(length(beh.time.trueTime), beh.nDLCpts*3);
beh.dlcAligned.tmp.aligned = zeros(length(behaviorOpts.timeTot), beh.nDLCpts*3);
beh.dlcAligned.data = zeros(length(behaviorOpts.timeTot)-behaviorOpts.bleachBuffer, beh.nDLCpts*3);
if beh.nDLCpts>0
    nms = beh.dlcData.Properties.VariableNames;
end

for j=1:beh.nDLCpts*3 % factor of 3 for x, y, and likelihood   
    dataChunk = beh.dlcData.(nms{j})(beh.time.imOn:beh.time.imOff);
    beh.dlcAligned.tmp.smooth(:,j) = smooth(dataChunk, beh.smoothing)';
    beh.dlcAligned.tmp.aligned(:,j) = interp1(beh.time.trueTime, beh.dlcAligned.tmp.smooth(:,j), behaviorOpts.timeTot,'nearest','extrap');    % interpolate to align time
    beh.dlcAligned.data(:,j) = beh.dlcAligned.tmp.aligned(behaviorOpts.bleachBuffer+1:end,j);
end
behAligned.dlcAligned.data = beh.dlcAligned.data;



function alignedBehaviorTot = concatenateBehaviorFiles(alignedBehaviorTot, alignedBehavior, time)

% same as standalone function concatenateBehaviorTrials.m, but without
% loading and without loop

if ~isfield(alignedBehaviorTot,'legVar')
    alignedBehaviorTot.legVar = alignedBehavior.legVar; %behSmooth; %
    alignedBehaviorTot.stim = alignedBehavior.stim;
    alignedBehaviorTot.drink = alignedBehavior.drink;
    alignedBehaviorTot.timeTot = time;
    alignedBehaviorTot.dlcData = alignedBehavior.dlcAligned.data;
else
    alignedBehaviorTot.legVar = [alignedBehaviorTot.legVar,alignedBehavior.legVar]; %behSmooth;
    alignedBehaviorTot.stim = [alignedBehaviorTot.stim,alignedBehavior.stim];
    alignedBehaviorTot.drink = [alignedBehaviorTot.drink,alignedBehavior.drink];
    alignedBehaviorTot.timeTot = [alignedBehaviorTot.timeTot,alignedBehaviorTot.timeTot(end)+time];
    alignedBehaviorTot.dlcData = [alignedBehaviorTot.dlcData; alignedBehavior.dlcAligned.data];
end


function parseStruct = getBehParsing(isImagingOn, traceFolder)
% converts a noisy measure of the indicator LED turning on/off in the
% behavior video into a binary vector by denoising, fitting a gaussian
% mixture model, then further smoothing the posterior to prevent flicker
imSm = TVL1denoise(double(isImagingOn), 0.1, 100);
if size(imSm,1)==1; imSm=imSm'; end
GMModel = fitgmdist(imSm,2);
P = posterior(GMModel,imSm);
[~,onIdx] = max(GMModel.mu);
[~,offIdx] = min(GMModel.mu);
gOn=TVL1denoise(P(:,onIdx), 0.1, 100); 
gOff=TVL1denoise(P(:,offIdx), 0.1, 100);
parseVec = gOn>gOff; % where posterior is higher for "on" distribution
dP = diff(parseVec);
starts = find(dP>0);
stops = find(dP<0);
parseStruct = struct('starts',starts,'stops',stops);
% if (min(diff(starts))<100) || (min(diff(stops))<100)
if (length(starts)>1) || (length(stops)>1)
    error('INDICATOR MEASURE IS FLICKERING. PARSING FAILURE');
elseif length(starts)~=length(stops)
    mkdir([traceFolder, 'crashReport']);
    crashMat = matfile([traceFolder, 'crashReport/','parseStruct.mat'],'Writable',true);
    crashMat.parseStruct = parseStruct;
    h = figure; 
    plot(isImagingOn);
    savefig(h, [traceFolder, 'crashReport/','isImagingOn.fig']);
    % save parseStruct and trace as fig
    error('NUMBER OF STARTS AND STOPS DO NOT MATCH. PARSING FAILURE');
end
