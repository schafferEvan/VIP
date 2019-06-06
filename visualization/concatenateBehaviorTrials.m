
function concatenateBehaviorTrials

addpath(genpath( '~/Dropbox/_code/'))
expDir = '/Users/evan/Dropbox/_sandbox/sourceExtraction/good/_runAndFeed/190422_f1/';
%expID = '0418'; %'0824_f3r1'; %'0824_f2r2_cur'; %'0828_f2r3'; %'0110'; %'0312_f4';

D = dir(expDir);
trialIndices = zeros(size(D));
for j=1:length(D)
    if D(j).isdir && strcmp(D(j).name(1),'r')
        trialIndices(j)=str2double(D(j).name(4));
    end
end
fileIdx = find(trialIndices);

alignedBehaviorTot = struct;

for j=1:length(fileIdx) % loop through trials to concatenate all data
    load([expDir,D(fileIdx(j)).name,'/alignedBehavAndStim.mat'])
    %load([expDir,'good/',expID,'/',D(fileIdx(j)).name,'/alignedBehavSmooth.mat'])

    if j==1
        alignedBehaviorTot.legVar = alignedBehavior.legVar; %behSmooth; %
        alignedBehaviorTot.stim = alignedBehavior.stim;
        alignedBehaviorTot.drink = alignedBehavior.drink;
        timeTot = time;
    else
        alignedBehaviorTot.legVar = [alignedBehaviorTot.legVar,alignedBehavior.legVar]; %behSmooth; %
        alignedBehaviorTot.stim = [alignedBehaviorTot.stim,alignedBehavior.stim];
        alignedBehaviorTot.drink = [alignedBehaviorTot.drink,alignedBehavior.drink];
        timeTot = [timeTot,timeTot(end)+time];
    end
end

alignedBehavior = alignedBehaviorTot;
time = timeTot;

save([expDir,'/alignedBehavAndStim.mat'],'alignedBehavior','time','-v7.3')