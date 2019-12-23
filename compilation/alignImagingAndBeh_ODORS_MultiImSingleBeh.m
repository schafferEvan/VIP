
function alignImagingAndBeh_ODORS_MultiImSingleBeh(codePath, traceFolder)%, repeatedBleachBuffer)


addpath(genpath(codePath))

repeatedBleachBuffer = 30;
imagingFile = [traceFolder,'F.mat'];
load(imagingFile,'trialFlag');

binDirectory = dir([traceFolder,'bin/fly*']); %dir([baseFolder,'scapeBehavior/',experiment,'*.mat']);

% get file order
tmp = zeros(size(binDirectory));
for j=1:length(binDirectory)
    tmp1 = strfind(binDirectory(j).name, 'run');
    tmp2 = strfind(binDirectory(j).name, '_');
    u = find(tmp2 > tmp1, 1);
    tmp(j)=str2double(binDirectory(j).name(tmp1+3:tmp2(u)-1));
end
if sum(isnan(tmp))>0; error('Invalid File Order'); end
[runIds,fileOrder] = sort(tmp,'ascend');


for j = 1:length(binDirectory)
    odorNum = runIds(j);
    binFilename = [binDirectory(fileOrder(j)).folder,'/',binDirectory(fileOrder(j)).name];
    
    
    fid = fopen(binFilename,'r');
    stimData = fread(fid, [8, inf], 'double');
    fclose(fid);
    PIDdata.time = stimData(1, :);
    PIDdata.data = stimData(6,:);
    PIDdata.odorNum = odorNum;

    
    
    infoDirectory = dir([traceFolder,'info/*.mat']);
    for k = 1:length(infoDirectory)
        if contains(infoDirectory(k).name,['run',num2str(odorNum)]) && (infoDirectory(k).name(1) ~= '.')
            infoFilename = [infoDirectory(k).folder,'/',infoDirectory(k).name];
            load(infoFilename,'info');
            break
        end
    end
    
    imRunLength = sum(trialFlag==odorNum);
    
    if j==1; bleachBuffer = 30*round(info.daq.scanRate);
    else;    bleachBuffer = repeatedBleachBuffer*round(info.daq.scanRate);
    end
    
    
    disp({'# scans: ',num2str(info.daq.numberOfScans);...
        'bleachBuffer: ',num2str(bleachBuffer);...
        'imRunLength: ',num2str(imRunLength)})
    if info.daq.numberOfScans-bleachBuffer~=imRunLength; error('length of traces does not match expected'); end
    
    
    timeTot = (1:info.daq.numberOfScans)/info.daq.scanRate;
    
    
    
    behaviorOpts.basicBehav = true; %legacy. ok if behav vid is const fps
    behaviorOpts.timeTot = timeTot;
    behaviorOpts.info = info;
    behaviorOpts.bleachBuffer = bleachBuffer;
    
    
    PIDAligned = formatPIDData(PIDdata, behaviorOpts);
    disp('PID file added')
    
    
    bmat = matfile([traceFolder,'alignedPID.mat'],'Writable',true);
    bmat.PIDAligned = PIDAligned;
end


function PIDAligned = formatPIDData(PIDdata, behaviorOpts)

PIDAligned.data = interp1(PIDdata.time, PIDdata.data, behaviorOpts.timeTot,'nearest','extrap');

PIDAligned.data = PIDAligned.data(behaviorOpts.bleachBuffer+1:end);
PIDAligned.time = behaviorOpts.timeTot;
PIDAligned.odorNum = PIDdata.odorNum;






