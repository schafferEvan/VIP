
function alignImagingAndBeh_LABEL_PROB_MultiImSingleBeh(codePath, traceFolder) %, repeatedBleachBuffer)

% ALIGNS *RAW* arHMM OUTPUT

% this is similar to alignImagingAndBehaviorMultiTrial.m, but whereas
% alignImagingAndBehaviorMultiTrial expects a directory of behavior files,
% this function expects a single behavior file containing multiple imaging
% runs.

% In 'experiment' folder, expects to find concatenated imaging file F.mat.
% Also expects a folder called 'behavior' containing each of the mat files
% ouput from extractBehaviorManual.m, named flyX_runY.mat, and a similarly
% populated folder called info

addpath(genpath(codePath))

repeatedBleachBuffer = 0; % this is now handled in postprocessing
% if nargin<3
%     repeatedBleachBuffer = 4; % default is that 4 sec of runs 2-end are ignored
% end
% if ischar(repeatedBleachBuffer); repeatedBleachBuffer = str2double(repeatedBleachBuffer); end % this can happen if called from bash script

%traceFolder = [baseFolder,experiment];%[baseFolder,'rustyOut/',experiment];
imagingFile = [traceFolder,'F.mat']; %[traceFolder,'all.mat'];
load(imagingFile,'trialFlag');
behaviorDirectory = dir([traceFolder,'behavior/fly*']); %dir([baseFolder,'scapeBehavior/',experiment,'*.mat']);
labelDirectory = dir([traceFolder,'beh_labels/*labels.mat']); %dir([baseFolder,'scapeBehavior/',experiment,'*.mat']);
%infoDirectory = dir([traceFolder,'info/fly*']); %dir([baseFolder,'scapeBehavior/',experiment,'*.mat']);
    
% get file order
% tmp = zeros(size(infoDirectory));
% for j=1:length(infoDirectory)
%     u = strfind(infoDirectory(j).name,'_');
%     r = strfind(infoDirectory(j).name,'run');
%     tmp(j)=str2double(infoDirectory(j).name(r+3:u(end)-1));
% end
[trials, ~, ~, runIds] = sortExperimentDirectory([traceFolder,'info/'],'_info');

if sum(isnan(runIds))>0; error('Invalid File Order'); end
% [runIds,fileOrder] = sort(tmp,'ascend');

labelsAligned = [];

behaviorFilename = [behaviorDirectory(1).folder,'/',behaviorDirectory(1).name];
behRaw = load(behaviorFilename);
labelFilename = [labelDirectory(1).folder,'/',labelDirectory(1).name];
labels = load(labelFilename);




% get list of imaging start and stop times from LED in behavior video
% if exists(crashreport) then parseStruct loads file, else runs function
crashname = dir([traceFolder, 'crashReport/']);
if ~isempty(crashname)
    crashFile = [traceFolder, 'crashReport/','parseStruct.mat'];
    load(crashFile,'parseStruct');
    behaviorOpts.parseStruct = parseStruct;
else
    parseStruct = getBehParsing(double(behRaw.traces.isImagingOn), traceFolder);
    if (length(parseStruct.starts)~=length(trials)) || (length(parseStruct.stops)~=length(trials))
        error('parse struct does not match info directory')
    end
    behaviorOpts.parseStruct = parseStruct;
end


% trials = trials(1:6);
for j=1:length(trials)
    %     if (j==4) || (j==5)
    %         continue
    %     end
    
    % load info file
    infoFile = [trials(j).folder,'/',trials(j).name];
    load(infoFile,'info');
   
    % grab length of correct piece of imaging data
    imRunLength = sum(trialFlag==runIds(j));
    
    % ignore first 30s of beginning of experiment, but not subsequent runs
    if j==1; bleachBuffer = 30*round(info.daq.scanRate);
    else;    bleachBuffer = repeatedBleachBuffer*round(info.daq.scanRate);
    end
    
    % check that extracted imaging data has correct number of frames
    disp({'# scans: ',num2str(info.daq.numberOfScans);...
        'bleachBuffer: ',num2str(bleachBuffer);...
        'imRunLength: ',num2str(imRunLength)})
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
    behaviorOpts.basicBehav = true; %legacy. ok if behav vid is const fps
    behaviorOpts.timeTot = timeTot;
    behaviorOpts.info = info;
    behaviorOpts.bleachBuffer = bleachBuffer;
    %time = timeTot(bleachBuffer+1:end);
    
    
    beh = formatLabelData(behRaw, labels, behaviorOpts, j);
    disp('Behavior file added')
    
    if j==1
        labelsAligned = beh.labels;
    else
        labelsAligned = [labelsAligned; beh.labels];
    end
end


bmat = matfile([traceFolder,'alignedLabels.mat'],'Writable',true);
bmat.labelsAligned = labelsAligned;




function behAligned = formatLabelData(beh, labels, behaviorOpts, iter)

% takes as input a matfile made by extractBehaviorManual.m or extractBehaviorAuto.m

beh.time.imOn = behaviorOpts.parseStruct.starts(iter);
beh.time.imOff = behaviorOpts.parseStruct.stops(iter);
% beh.smoothing = max(round( ((beh.time.imOff-beh.time.imOn)/behaviorOpts.info.daq.scanLength)/behaviorOpts.info.daq.scanRate ),1);    % smooth by behavior frames per imaging frame
beh.time.trueTime = (1:(beh.time.imOff-beh.time.imOn+1))/(beh.time.imOff-beh.time.imOn+1)*behaviorOpts.info.daq.scanLength;            % fixed?

% align labels using windowed mean of original
labelsWindowed = labels.states(beh.time.imOn:beh.time.imOff,:);

NB = size(labels.states,2);
beh.traces.fullLabels = zeros(length(behaviorOpts.timeTot), NB);
[~,strt] = find(beh.time.trueTime<.5*behaviorOpts.timeTot(1),1,'last');
[~,stp] = find(beh.time.trueTime<.5*sum(behaviorOpts.timeTot(1:2)),1,'last');
beh.traces.fullLabels(1,:) = mean(labelsWindowed(strt:stp,:),1);

for j=2:length(behaviorOpts.timeTot)-1
    strt = stp+1;
    [~,stp] = find(beh.time.trueTime<.5*sum(behaviorOpts.timeTot(j:j+1)),1,'last');
    beh.traces.fullLabels(j,:) = mean(labelsWindowed(strt:stp,:),1);
end

strt = stp+1;
[~,stp] = find(beh.time.trueTime<behaviorOpts.timeTot(end),1,'last');
beh.traces.fullLabels(j,:) = mean(labelsWindowed(strt:stp,:),1);

behAligned.labels = beh.traces.fullLabels(behaviorOpts.bleachBuffer+1:end,:);






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
if (min(diff(starts))<100) || (min(diff(stops))<100)
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


