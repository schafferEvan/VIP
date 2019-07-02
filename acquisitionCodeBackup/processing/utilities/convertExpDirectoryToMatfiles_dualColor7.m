
% reads a directory of folders, each containing tiff stacks, and exports a
% set of matfiles containing chunks of data and a template


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flyNum = '1';
%flyNum = 2;
chunkSize = 300;%201; %871;%677;%751;                                                    % number of time steps per matfile % make it bigger (around 300)
zeroThreshold = 0; %120;                                                % threshold for re-zeroing data % don't change
cropping.x = 28:154; %28:110;%31:121;%26:118; %28:119; %32:129;%32:129;%41:124;%35:117;%23:104;%28:114;%24:117;%35:127;%24:113;%45:124;%25:117;%38:120;%36:128;%39:132;%32:164;%25:180;%40:195; %25:145; %43:115;%36:119;%39:167;
cropping.y = 62:327; %73:336;%67:335;%57:323; %78:351; %42:321;%42:321;%15:298;%57:309;%77:362;%36:319; %50:322;%53:357;%44:357;%32:326;%54:354;%60:330;%39:304;%63:346;%58:339;%25:360;%75:350; %60:333; %590:862;%12:290;%[];
cropping.z = 8:185; %1:165;%5:189;%10:139; %19:158; %10:186; %10:186; %17:152;%22:159;%12:147;%18:152;%7:144;%37:182;%12:127;%20:216;%16:141;%31:216;%10:186;%18:205;%42:246;%1:204;%12:190; %7:180; %13:149;%42:184;%29:149;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ismac
    addpath(genpath('/Volumes/dropboxdrive 1/Dropbox/05_GitHub/NoRMCorre-master/NoRMCorre-master-new/'))
    addpath(genpath('/Volumes/SCAPEdata1/code/pipeline/'))
    baseFolderName = '/Volumes/workdrive/';
    experimentName = '2019_05_03_57C10_NLS6f_Dilp5_GtACR_noAtRctrl_walk';%'2018_07_19_hungry';
    datafolder = [baseFolderName,experimentName,'/tiff_stacks/'];     % Folder containing folder of tiff stacks
    savePath = ['/Volumes/workdrive3/',experimentName,'/_matfiles/'];
    infoFile = dir([datafolder,'fly',flyNum,'_info/*.mat']);%'/Volumes/dataFast/habaRegistered/2018_08_24_odor/mats/new/fly3run2/'; %'/Users/evan/Desktop/hungerRaw/';%'/Volumes/SCAPEdata1/scratchData/2018_08_01_IRtest/matfiles/registered/';%'/Volumes/data/_scape/data/_outMats/'; %
    %infoFile = dir([datafolder,'fly',num2str(flyNum),'_info/*.mat']);
    if length(infoFile)>1; infoFile=infoFile(end); end
    mkdir(savePath);
else
    addpath(genpath('C:\Users\Axel-SCAPE\Desktop\motion_correction\motion_correction'))
    baseFolderName = 'H:\';%'D:\';
    experimentName = '2019_06_30_Nsyb_NLS6s_walk';%'2019_05_16_57C10_NLS6f_sNPF_red_Su_walk';%'2018_07_19_hungry';
    expFolder = [baseFolderName,experimentName,'\'];
    datafolder = [baseFolderName,experimentName,'\tiff_stacks\'];     % Folder containing folder of tiff stacks
    savePath = [datafolder,'_matfiles_new\'];
    %infoFile = dir([datafolder,'fly',num2str(flyNum),'_info\*.mat']);%'/Volumes/dataFast/habaRegistered/2018_08_24_odor/mats/new/fly3run2/'; %'/Users/evan/Desktop/hungerRaw/';%'/Volumes/SCAPEdata1/scratchData/2018_08_01_IRtest/matfiles/registered/';%'/Volumes/data/_scape/data/_outMats/'; %
    infoFile = dir([expFolder,'fly',flyNum,'*_info.mat']);%dir([datafolder,'fly',flyNum,'_info\*.mat']);
    if length(infoFile)>1; infoFile=infoFile(end); end
    mkdir(savePath);
end
load([infoFile.folder,'\',infoFile.name]);
timeToExclude = 30*round(info.daq.scanRate);                                              % number of time steps at beginning to skip %framerate * 30 seconds


%d = dir([datafolder,'fly',flyNum,'*']);    % list of datasets (folders)
d = dir([datafolder,'fly',num2str(flyNum),'_run*']);
%d = dir([datafolder,'fly',num2str(flyNum),'*']);
templateID = ones(size(d));                                           % specify where to get the template. set as 1:length(d) to give each experiment its own template.  Set as ones(size(d)) to use a single template

% m = matfile('E:\SCAPE_DATA\2018_08_16_IRtest\tiff_stacks\_matfiles\2018_08_16_IRtest_fly1_run1_run1_301.mat');
% template = m.template;

for j=1:length(d)
    trialfolder = [datafolder,d(j).name,'\'];
    template = [];
%     if j == 1
%         makeManyMatfilesFromTiffStack_dualColor(trialfolder,savePath,experimentName,chunkSize,zeroThreshold,timeToExclude,cropping);
%     else
        makeManyMatfilesFromTiffStack_dualColor(trialfolder,savePath,experimentName,chunkSize,zeroThreshold,0,cropping,template);
%     end
end