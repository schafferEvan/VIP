
function retroactivelyIgnoreTrials(traceFolder, trimTrials)
% removes trials/runs from further analysis (useful if something goes wrong
% late in an experiment)
% trimTrials is a vector of trial *names*, i.e. [3, 4.1, 5]

% NOTE: WHEREAS RETROACTIVELYADBLEACHBUFFER.M IS MEANT TO BE RUN BEFORE
% ALIGNING IMAGING AND BEHAVIOR DATA, THIS SCRIPT IS MEANT TO RUN AFTER
% (BUT STILL BEFORE POSTPROCESSING.PY)

%if ischar(trimSecs); trimSecs = str2double(trimSecs); end % this can happen if called from bash script
disp(['Removing runs ',num2str(trimTrials)])

addpath(genpath('../..'))

mkdir([traceFolder,'untrimmedTraces/allRuns/'])


% -------- from green ----------------
file = 'F.mat';
disp(['Trimming ',file])

copyfile([traceFolder,file],[traceFolder,'untrimmedTraces/allRuns/',file])
load([traceFolder,file]);
[~, ~, ~, runIds] = sortExperimentDirectory([traceFolder,'info/'],'_info');
%load([traceFolder,'info/',trials(1).name]);

trUnique = unique(trialFlag);
trialFlagNew = trialFlag;

for j=1:length(runIds)
    if trUnique(j)~=runIds(j); error( 'trial flags do not match' ); end
    if isempty(find(trimTrials==runIds(j),1))
        continue; % continue if this run isn't being removed
    end
    thisRun = (trialFlag==runIds(j));
    trialFlagNew(thisRun) = NaN;
end
toKeep = find(~isnan(trialFlagNew));

F = F(:,toKeep);
FR = FR(:,toKeep);
trialFlag = trialFlagNew(~isnan(trialFlagNew));
save([traceFolder,file],'A','F','FR','traceFolder','trialFlag' );



% -------- from red ----------------
file = 'F_fromRed.mat';
disp(['Trimming ',file])

copyfile([traceFolder,file],[traceFolder,'untrimmedTraces/allRuns/',file])
load([traceFolder,file]);

F = F(:,toKeep);
FR = FR(:,toKeep);
trialFlag = trialFlagNew(~isnan(trialFlagNew));
save([traceFolder,file],'A','F','FR','traceFolder','trialFlag' );




% -------- aligned behavior ----------------
file = 'alignedBehavAndStim.mat';
disp(['Trimming ',file])

copyfile([traceFolder,file],[traceFolder,'untrimmedTraces/allRuns/',file])
load([traceFolder,file]);

time = time(:,toKeep);
alignedBehavior.legVar = alignedBehavior.legVar(toKeep);
alignedBehavior.stim = alignedBehavior.stim(toKeep);
alignedBehavior.drink = alignedBehavior.drink(toKeep);
alignedBehavior.timeTot = alignedBehavior.timeTot(toKeep);
alignedBehavior.dlcData = alignedBehavior.dlcData(toKeep,:);
save([traceFolder,file],'alignedBehavior','time');

