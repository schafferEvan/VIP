
function retroactivelyAddBleachBuffer(traceFolder, trimSecs)
% allows for a bleach buffer to be added to F, FR, and behavior 
% after the source extraction steps have been run.

if ischar(trimSecs); trimSecs = str2double(trimSecs); end % this can happen if called from bash script
disp(['Trimming an additional ',num2str(trimSecs), 'sec for bleaching'])

addpath(genpath('../..'))

if isfolder([traceFolder,'untrimmedTraces/'])
    if exist([traceFolder,'untrimmedTraces/F.mat'],'file')
        disp(' Trimming already happened.  Skipping.')
        return
    end
end
mkdir([traceFolder,'untrimmedTraces/'])


% -------- from green ----------------
file = 'F.mat';
disp(['Trimming ',file])

copyfile([traceFolder,file],[traceFolder,'untrimmedTraces/',file])
load([traceFolder,file]);
[trials, ~, ~, runIds] = sortExperimentDirectory([traceFolder,'info/'],'_info');
load([traceFolder,'info/',trials(1).name]);

bleachBufferNew = trimSecs*round(info.daq.scanRate);

trUnique = unique(trialFlag);
trialFlagNew = trialFlag;

for j=1:length(runIds)
    if trUnique(j)~=runIds(j); error( 'trial flags do not match' ); end
    if j==1
        continue; %don't do anything to first run
    end
    thisRun = find(trialFlag==runIds(j));
    trialFlagNew(thisRun(1:bleachBufferNew)) = NaN;
end
toKeep = find(~isnan(trialFlagNew));

F = F(:,toKeep);
FR = FR(:,toKeep);
trialFlag = trialFlagNew(~isnan(trialFlagNew));
save([traceFolder,file],'A','F','FR','traceFolder','trialFlag' );



% -------- from red ----------------
file = 'F_fromRed.mat';
disp(['Trimming ',file])

copyfile([traceFolder,file],[traceFolder,'untrimmedTraces/',file])
load([traceFolder,file]);

F = F(:,toKeep);
FR = FR(:,toKeep);
trialFlag = trialFlagNew(~isnan(trialFlagNew));
save([traceFolder,file],'A','F','FR','traceFolder','trialFlag' );

