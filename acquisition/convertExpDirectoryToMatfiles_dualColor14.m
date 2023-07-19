
% reads a directory of folders, each containing tiff stacks, and exports a
% set of matfiles containing chunks of data and a template


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flyNum = '2';
%flyNum = 2;
chunkSize = 200;%201; %871;%677;%751;                                                    % number of time steps per matfile % make it bigger (around 300)
zeroThreshold = 0; %120;                                                % threshold for re-zeroing data % don't change
cropping.x = 50:130;%33:123;%43:142;%32:131;%43:147;%37:122;%25:118;%41:116;%25:121;%32:120;%38:118;%34:128;%38:116;%23:114;%33:128;%40:131;%41:128;%37:114;%25:113;%39:117;%26:121;%40:126;%33:129;%35:130;%43:131;%33:120;%37:113;%30:106;%28:108; %26:118; %28:119; %32:129;%32:129;%41:124;%35:117;%23:104;%28:114;%24:117;%35:127;%24:113;%45:124;%25:117;%38:120;%36:128;%39:132;%32:164;%25:180;%40:195; %25:145; %43:115;%36:119;%39:167;
cropping.y = 85:225;%91:348; %68:360;%44:376;%82:350;%64:334;%54:357;%65:348;%70:318;%56:326;%70:312;%57:327;%37:296;%76:338;%61:327;%64:325; %73:338;%67:316;%58:329;%10:299;%62:353; %41:360;%78:339;%60:310;%77:341;%51:336;%50:334;%77:339;%82:349; %57:323; %78:351; %42:321;%42:321;%15:298;%57:309;%77:362;%36:319; %50:322;%53:357;%44:357;%32:326;%54:354;%60:330;%39:304;%63:346;%58:339;%25:360;%75:350; %60:333; %590:862;%12:290;%[];
cropping.z = 80:160;%30:140;%27:196;%13:150;%35:146;%22:171;%1:138;%7:160;%8:166;%15:195;%14:139;%18:160;%22:140;%18:187;%11:193;%25:191;%25:175;%39:198;%17:191;%14:139;%8:185;%35:172; %18:170;%25:200;%39:163;%25:217;%13:157;%32:166;%45:176;%23:138; %10:139; %19:158; %10:186; %10:186; %17:152;%22:159;%12:147;%18:152;%7:144;%37:182;%12:127;%20:216;%16:141;%31:216;%10:186;%18:205;%42:246;%1:204;%12:190; %7:180; %13:149;%42:184;%29:149;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ismac
    addpath(genpath('/Users/evan/Dropbox/_code/'))
    baseFolderName = '/Volumes/RASPIUSB/_scapeTemp/';
    experimentName = '2019_05_14_57C10_NLS6f_CrzDh44_red_Su_walk';%'2018_07_19_hungry';
    expFolder = [baseFolderName,experimentName,'/'];
    datafolder = [expFolder,'/tiff_stacks/'];     % Folder containing folder of tiff stacks
    savePath = [baseFolderName,experimentName,'/_matfiles/'];
    infoFile = dir([expFolder,'fly',flyNum,'*_info.mat']);%'/Volumes/dataFast/habaRegistered/2018_08_24_odor/mats/new/fly3run2/'; %'/Users/evan/Desktop/hungerRaw/';%'/Volumes/SCAPEdata1/scratchData/2018_08_01_IRtest/matfiles/registered/';%'/Volumes/data/_scape/data/_outMats/'; %
    %infoFile = dir([datafolder,'fly',num2str(flyNum),'_info/*.mat']);
    if length(infoFile)>1; infoFile=infoFile(end); end
    mkdir(savePath);
    load([infoFile.folder,'/',infoFile.name]);
else
    addpath(genpath('C:\Users\Axel-SCAPE\Desktop\motion_correction\motion_correction'))
    addpath(genpath('C:\Users\Axel-SCAPE\Documents\_code\analysis\utilities'))
    baseFolderName = 'D:\';%'D:\';
    experimentName = '2022_01_08_oviDN_NLS6s';%'2019_05_16_57C10_NLS6f_sNPF_red_Su_walk';%'2018_07_19_hungry';
    expFolder = [baseFolderName,experimentName,'\'];
    datafolder = [baseFolderName,experimentName,'\tiff_stacks\'];     % Folder containing folder of tiff stacks
    savePath = [datafolder,'_matfiles_new\'];
    infoFile = dir([expFolder,'fly',flyNum,'*_info.mat']);%dir([datafolder,'\fly',num2str(flyNum),'_info\*_info.mat']);%'/Volumes/dataFast/habaRegistered/2018_08_24_odor/mats/new/fly3run2/'; %'/Users/evan/Desktop/hungerRaw/';%'/Volumes/SCAPEdata1/scratchData/2018_08_01_IRtest/matfiles/registered/';%'/Volumes/data/_scape/data/_outMats/'; %
    %infoFile = dir(c);
    if length(infoFile)>1; infoFile=infoFile(end); end
    mkdir(savePath);
    load([infoFile.folder,'\',infoFile.name]);
end

timeToExclude = 30*round(info.daq.scanRate);                                              % number of time steps at beginning to skip %framerate * 30 seconds


%d = dir([datafolder,'fly',flyNum,'*']);    % list of datasets (folders)
d = dir([datafolder,'fly',num2str(flyNum),'_run*']);
%d = dir([datafolder,'fly',num2str(flyNum),'*']);
templateID = ones(size(d));                                           % specify where to get the template. set as 1:length(d) to give each experiment its own template.  Set as ones(size(d)) to use a single template

% m = matfile('E:\SCAPE_DATA\2018_08_16_IRtest\tiff_stacks\_matfiles\2018_08_16_IRtest_fly1_run1_run1_301.mat');
% template = m.template;

for j=1:length(d)
    if ismac
        trialfolder = [datafolder,d(j).name,'/'];
    else
        trialfolder = [datafolder,d(j).name,'\'];
    end
    template = [];
    if j == 1
        makeManyMatfilesFromTiffStack_dualColor(trialfolder,savePath,experimentName,chunkSize,zeroThreshold,timeToExclude,cropping);
    else
        makeManyMatfilesFromTiffStack_dualColor(trialfolder,savePath,experimentName,chunkSize,zeroThreshold,0,cropping,template);
    end
end