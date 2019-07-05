function extract_F_from_conComp(codePath, experimentFolder)
% 
% % from connected components structure generated by watershed filter on
% % static image, reads directory of original matfiles and extracts raw F
% 
% %clear
% %pause(60*100)
% 
% if ismac
%     addpath(genpath('~/Dropbox/_AxelLab/matlab/calcium-signal-extraction/'))
%     addpath(genpath( '~/Dropbox/_code/pipeline/'))
%     experimentFolder = '/Volumes/SCAPEdata1/finalData/2019_03_12_Nsyb_NLS6s_Su/fly4/'; %2018_11_01_looming/sample/2018_11_01_looming/';%'/Volumes/SCAPEdata1/scratchData/2018_10_31_IRtest/_matfiles_pt1/';%1024mcfo/';%'/Volumes/data/registered/ir0831/matfiles/';%'/Volumes/dataFast/habaRegistered/2018_08_24_odor/mats/new/fly3run2/'; %'/Users/evan/Desktop/hungerRaw/';%'/Volumes/SCAPEdata1/scratchData/2018_08_01_IRtest/matfiles/registered/';%'/Volumes/data/_scape/data/_outMats/'; %
%     %infoFile = dir([experimentFolder,'info/*.mat']);%'/Volumes/dataFast/habaRegistered/2018_08_24_odor/mats/new/fly3run2/'; %'/Users/evan/Desktop/hungerRaw/';%'/Volumes/SCAPEdata1/scratchData/2018_08_01_IRtest/matfiles/registered/';%'/Volumes/data/_scape/data/_outMats/'; %
%     %if length(infoFile)>1; infoFile=infoFile(end); end
%     savePath = '~/Dropbox/_AxelLab/_data/_scape/movies/';
% elseif isunix
%     addpath(genpath('/home/analysis-pc/00_Analysis/motion_correction/'))
%     addpath(genpath('/home/analysis-pc/00_Analysis/calcium-signal-extraction/'))
%     experimentFolder = '/home/analysis-pc/rawData/20171025_nSyb_fly2_reg/registered/reregistered/';
%     savePath = '/home/analysis-pc/rawData/movies/';
% else
%     addpath(genpath('E:\Dropbox\GitHub\calcium-signal-extraction'))
%     addpath(genpath('E:\Dropbox\GitHub\eftyMotionCorr'))
%     experimentFolder = 'D:\SCAPEdataworkingfolder\_outMats\';
%     savePath = 'D:\SCAPEdataworkingfolder\movies\';
% end
% 
% load([experimentFolder,'Yproj/cc.mat'],'cc')
% 
% trials = dir([experimentFolder,'*.mat']); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%dir([experimentFolder,'*.h5']);
% frameNum = zeros(size(trials));
% runNum = zeros(size(trials));
% for j=1:length(trials)
%     runLoc = strfind(trials(j).name,'run');
%     underscoreLoc = strfind(trials(j).name,'_');
%     regLoc = strfind(trials(j).name,'reg');
%     if isempty(regLoc)
%         regLoc = strfind(trials(j).name,'.');
%     end
%     if ~isempty(runLoc)
%         uRel = find(underscoreLoc>runLoc,1,'first');
%         runNum(j) = str2double(trials(j).name(runLoc+3:underscoreLoc(uRel)-1));
%     else
%         runNum(j) = 1;
%     end
%     frameNum(j) = str2double(trials(j).name(underscoreLoc(end)+1:regLoc(1)-1));     
% end
% [~,trialOrder] = sort( runNum*10^ceil(log10(max(frameNum))) + frameNum, 'ascend' );
% 

addpath(genpath(codePath))
[trials, trialOrder, runNum, runNumConvert] = sortExperimentDirectory(experimentFolder);

save([experimentFolder,'Yproj/Ysum.mat'],'Ysum','Rsum')
save([experimentFolder,'Yproj/cc.mat'],'Ycc','Rcc')

%% extract F for all cells, all files
trialPath = [experimentFolder,trials(trialOrder(1)).name];
m = matfile(trialPath);
tStepsReg = size(m,'Y',4);
trialPath = [experimentFolder,trials(trialOrder(end)).name];
m = matfile(trialPath);
%tStepsEnd = size(m,'Y',4);

% get last piece of each trial
runListTmp = unique(runNumConvert);
tStepsEnd = 0;
for j=1:length(runListTmp)
    tmp = find(runNumConvert==runListTmp(j),1,'last');
    trialPath = [experimentFolder,trials(trialOrder(tmp)).name];
    m = matfile(trialPath);
    tStepsEnd = tStepsEnd + size(m,'Y',4);
end

% parsing from green -----------------------------------------------------
A = sparse(prod(Ycc.ImageSize), Ycc.NumObjects);
for j=1:Ycc.NumObjects
    disp(j)
    A(Ycc.PixelIdxList{j},j)=1;
end

F = zeros(Ycc.NumObjects, (length(trials)-length(runListTmp))*tStepsReg + tStepsEnd );
FR = zeros(Ycc.NumObjects, (length(trials)-length(runListTmp))*tStepsReg + tStepsEnd );

% parsing from red -------------------------------------------------------

if ~isempty(Rcc)
    Ar = sparse(prod(Rcc.ImageSize), Rcc.NumObjects);
    for j=1:Rcc.NumObjects
        disp(j)
        Ar(Rcc.PixelIdxList{j},j)=1;
    end
    Fr = zeros(Rcc.NumObjects, (length(trials)-length(runListTmp))*tStepsReg + tStepsEnd );
    FRr = zeros(Rcc.NumObjects, (length(trials)-length(runListTmp))*tStepsReg + tStepsEnd );
else
    Ar = [];
    Fr = [];
    FRr = [];
end

trialFlag = zeros((length(trials)-length(runListTmp))*tStepsReg + tStepsEnd, 1);

% ------------------------------------------------------------------------
t=0;

for ii=1:length(trials)
    
    trialPath = [experimentFolder,trials(trialOrder(ii)).name];
    display(['file ',num2str(ii),': ',trials(trialOrder(ii)).name])
    m = matfile(trialPath);
    Y = m.Y;
    
    T = size(Y,4);
    Y = reshape(Y,prod(Ycc.ImageSize),T);
    F(:,t+(1:T)) = A'*Y;
    if ~isempty(Ar)
        Fr(:,t+(1:T)) = Ar'*Y;
    end
    trialFlag(t+(1:T))=runNumConvert(trialOrder(ii));
    clear Y
    
    if ~isempty(Ar)
        R = m.R;
        R = reshape(R,prod(Ycc.ImageSize),T);
        FR(:,t+(1:T)) = A'*R;
        FRr(:,t+(1:T)) = Ar'*R;
    end
    clear R
    
    t = t+T;
    
end

save([experimentFolder,'Yproj/F.mat'],'F','A','FR','trialFlag')

A = Ar;
F = Fr;
FR = FRr;
save([experimentFolder,'Yproj/F_fromRed.mat'],'F','A','FR','trialFlag')
