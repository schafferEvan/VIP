

function compile_sumImage(codePath, experimentFolder)

addpath(genpath(codePath))
    
% clear
% %pause(60*150)
% 
% if ismac
%     addpath(genpath('~/Dropbox/_AxelLab/matlab/calcium-signal-extraction/'))
%     addpath(genpath( '~/Dropbox/_code/pipeline/'))
%     experimentFolder = '/Volumes/SCAPEdata1/finalData/2019_05_01_57C10_NLS6f_Dilp5_GtACR/';%'/Volumes/SCAPEdata1/finalData/2019_05_15_57C10_NLS6f_Oct_red_Su_walk/';%'/Volumes/SCAPEdata1/finalData/2019_04_22_Nsyb_NLS6s_Su_walk/fly1/run4/'; %'/Volumes/dataFast/habaRegistered/2019_02_14_Nsyb_NLS6s_Su/';%'/Volumes/dataFast/habaRegistered/2018_11_20_IRtest/';%'/Volumes/dataFast2/habaRegistered/2018_08_24_odor/fly3run2/'; %2018_11_01_looming/sample/2018_11_01_looming/';%'/Volumes/SCAPEdata1/scratchData/2018_10_31_IRtest/_matfiles_pt1/';%1024mcfo/';%'/Volumes/data/registered/ir0831/matfiles/';%'/Volumes/dataFast/habaRegistered/2018_08_24_odor/mats/new/fly3run2/'; %'/Users/evan/Desktop/hungerRaw/';%'/Volumes/SCAPEdata1/scratchData/2018_08_01_IRtest/matfiles/registered/';%'/Volumes/data/_scape/data/_outMats/'; %
%     %infoFile = dir([experimentFolder,'info/*.mat']);%'/Volumes/dataFast/habaRegistered/2018_08_24_odor/mats/new/fly3run2/'; %'/Users/evan/Desktop/hungerRaw/';%'/Volumes/SCAPEdata1/scratchData/2018_08_01_IRtest/matfiles/registered/';%'/Volumes/data/_scape/data/_outMats/'; %
%     %if length(infoFile)>1; infoFile=infoFile(end); end
%     %savePath = '~/Dropbox/_AxelLab/_data/_scape/movies/';
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

[trials, trialOrder] = sortExperimentDirectory(experimentFolder);

% more extensive check of whether files are corrupt
for ii=1:length(trials)    
    trialPath = [experimentFolder,trials(trialOrder(ii)).name];

    try
        m = matfile(trialPath);
        sy = size(m,'Y');
        sr = size(m,'R');
        m.Y(sy(1),sy(2),sy(3),sy(4));
        m.R(sr(1),sr(2),sr(3),sr(4));
    catch 
        disp([trials(trialOrder(ii)).name,' is corrupted']);
    end
end

%% sum image calc

for ii=1:length(trials)
    
    trialPath = [experimentFolder,trials(trialOrder(ii)).name];
    display(['file ',num2str(ii),': ',trials(trialOrder(ii)).name])
    %m = matfile(trialPath);
    %Y = m.Y;
    try
        load(trialPath);
    catch % rare issue with partially readable matfile
        m = matfile(trialPath);
        Y = m.Y;
        R = m.R;
    end
    if ii==1
        Ysum = max(Y,[],4);
        if exist('R','var')
            Rsum = max(R,[],4);
        else
            Rsum = [];
        end
    else
        Ysum = max(Ysum, max(Y,[],4));
        if exist('R','var')
            Rsum = max(Rsum, max(R,[],4));
        else
            Rsum = [];
        end
    end
    
end

mkdir([experimentFolder,'Yproj'])
save([experimentFolder,'Yproj/Ysum.mat'],'Ysum','Rsum')

[im_bw_out,Ycc,regionProps] = segmentHessian_SCAPE_new(Ysum);
if ~isempty(Rsum)
    [~,Rcc,~] = segmentHessian_SCAPE_new(Rsum);
else
    Rcc = [];
end
save([experimentFolder,'Yproj/cc.mat'],'Ycc','Rcc')