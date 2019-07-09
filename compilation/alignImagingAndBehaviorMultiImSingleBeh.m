
function alignImagingAndBehaviorMultiImSingleBeh(codePath, traceFolder)
% this is similar to alignImagingAndBehaviorMultiTrial.m, but whereas
% alignImagingAndBehaviorMultiTrial expects a directory of behavior files,
% this function expects a single behavior file containing multiple imaging
% runs.

% In 'experiment' folder, expects to find concatenated imaging file F.mat.
% Also expects a folder called 'behavior' containing each of the mat files
% ouput from extractBehaviorManual.m, named flyX_runY.mat, and a similarly
% populated folder called info

addpath(genpath(codePath))

%traceFolder = [baseFolder,experiment];%[baseFolder,'rustyOut/',experiment];
imagingFile = [traceFolder,'F.mat']; %[traceFolder,'all.mat'];
load(imagingFile,'trialFlag');
behaviorDirectory = dir([traceFolder,'behavior/fly*']); %dir([baseFolder,'scapeBehavior/',experiment,'*.mat']);
infoDirectory = dir([traceFolder,'info/fly*']); %dir([baseFolder,'scapeBehavior/',experiment,'*.mat']);
    
% get file order
% tmp = zeros(size(infoDirectory));
% for j=1:length(infoDirectory)
%     u = strfind(infoDirectory(j).name,'_');
%     r = strfind(infoDirectory(j).name,'run');
%     tmp(j)=str2double(infoDirectory(j).name(r+3:u(end)-1));
% end
[~, fileOrder, ~, runIds] = sortExperimentDirectory([traceFolder,'info/'],'_info');
if sum(isnan(runIds))>0; error('Invalid File Order'); end
% [runIds,fileOrder] = sort(tmp,'ascend');

alignedBehaviorTot = struct;

behaviorFilename = [behaviorDirectory(1).folder,'/',behaviorDirectory(1).name];
behRaw = load(behaviorFilename);

% get list of imaging start and stop times from LED in behavior video
parseStruct = getBehParsing(double(behRaw.traces.isImagingOn));
if (length(parseStruct.starts)~=length(infoDirectory)) || (length(parseStruct.stops)~=length(infoDirectory))
    error('parse struct does not match info directory')
end
behaviorOpts.parseStruct = parseStruct;


for j=1:length(infoDirectory)
    
    index = fileOrder(j);
    % load info file
    infoFile = [infoDirectory(index).folder,'/',infoDirectory(index).name];
    load(infoFile,'info');
   
    % grab length of correct piece of imaging data
    imRunLength = sum(trialFlag==runIds(index));
    
    % ignore first 30s of beginning of experiment, but not subsequent runs
    if j==1; bleachBuffer = 30*round(info.daq.scanRate);
    else;    bleachBuffer = 0;
    end
    
    % check that extracted imaging data has correct number of frames
    if info.daq.numberOfScans-bleachBuffer~=imRunLength; error('length of traces does not match expected'); end
    
    % populate imaging time vector
    timeTot = (1:info.daq.numberOfScans)/info.daq.scanRate;
    timeTmp = timeTot(bleachBuffer+1:end);
    %     if j==1
    %         timeTot = timeTot; %1/info.daq.scanRate:1/info.daq.scanRate:info.daq.scanLength;
    %     else
    %         timeTot = [timeTot, timeTot(end) + timeTot]; %1/info.daq.scanRate:1/info.daq.scanRate:info.daq.scanLength;
    %     end
    
    
    
    % parse behavior
    behaviorOpts.basicBehav = true; %false; % true => old defaults, ok if behav is const 30fps
    behaviorOpts.timeTot = timeTot;
    behaviorOpts.info = info;
    behaviorOpts.bleachBuffer = bleachBuffer;
    %time = timeTot(bleachBuffer+1:end);
    
    
    beh = formatBehaviorData(behRaw, behaviorOpts, index);
    disp(['Behavior file added: ', index]);
    alignedBehaviorTot = concatenateBehaviorFiles(alignedBehaviorTot, beh, timeTmp); % ********** this needs to change

end

    

bmat = matfile([traceFolder,'alignedBehavAndStim.mat'],'Writable',true);
bmat.alignedBehavior = alignedBehaviorTot;
bmat.time = alignedBehaviorTot.timeTot;




function behAligned = formatBehaviorData(beh, behaviorOpts, iter)

% takes as input a matfile made by extractBehaviorManual.m or extractBehaviorAuto.m

beh.time.imOn = behaviorOpts.parseStruct.starts(iter);
beh.time.imOff = behaviorOpts.parseStruct.stops(iter);
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



function parseStruct = getBehParsing(isImagingOn)
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
if (min(diff(starts))<100) || (min(diff(stops))<100)
    error('INDICATOR MEASURE IS FLICKERING. PARSING FAILURE');
elseif length(starts)~=length(stops)
    error('NUMBER OF STARTS AND STOPS DO NOT MATCH. PARSING FAILURE');
end
parseStruct = struct('starts',starts,'stops',stops);

