
function alignImagingAndBehaviorMultiImSingleBeh(codePath, traceFolder) %, repeatedBleachBuffer)
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

alignedBehaviorTot = struct;

behaviorFilename = [behaviorDirectory(1).folder,'/',behaviorDirectory(1).name];
behRaw = load(behaviorFilename);

% check for deeplabcut output and load if available
dlcname = dir([traceFolder,'behavior/*DeepCut*.csv']);
if ~isempty(dlcname)
    behRaw.nDLCpts = 8; disp('update dlcRead to allow for variable number of points')
    behRaw.dlcData = dlcRead([traceFolder, 'behavior/', dlcname.name],behRaw.nDLCpts);
else
    behRaw.nDLCpts = 0; 
    behRaw.dlcData = [];
end

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
    if j==1
        bleachBuffer = 30*round(info.daq.scanRate); % 30*info.daq.scanRate; % use float version for true time
        %bleachBufferTrunc = 30*round(info.daq.scanRate); % use approx version for sync check
    else    
        bleachBuffer = repeatedBleachBuffer*round(info.daq.scanRate);
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
    
    
    beh = formatBehaviorData(behRaw, behaviorOpts, j);
    disp('Behavior file added')
    alignedBehaviorTot = concatenateBehaviorFiles(alignedBehaviorTot, beh, timeTmp);

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
if ~isempty(beh.dlcData)
    nms = beh.dlcData.Properties.VariableNames;
    
    for j=1:beh.nDLCpts*3 % factor of 3 for x, y, and likelihood
        dataChunk = beh.dlcData.(nms{j})(beh.time.imOn:beh.time.imOff);
        beh.dlcAligned.tmp.smooth(:,j) = smooth(dataChunk, beh.smoothing)';
        beh.dlcAligned.tmp.aligned(:,j) = interp1(beh.time.trueTime, beh.dlcAligned.tmp.smooth(:,j), behaviorOpts.timeTot,'nearest','extrap');    % interpolate to align time
        beh.dlcAligned.data(:,j) = beh.dlcAligned.tmp.aligned(behaviorOpts.bleachBuffer+1:end,j);
    end
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
    alignedBehaviorTot.legVar = [alignedBehaviorTot.legVar,alignedBehavior.legVar]; %behSmooth; %
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
if (length(starts)>1) || (length(stops)>1)
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
end


