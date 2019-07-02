
% reads a directory of folders, each containing tiff stacks, and exports a
% set of matfiles containing chunks of data and a template


%pause(4*60*60);


if ismac
    addpath(genpath('/Users/evan/Dropbox/_code/NoRMCorre'))
    baseFolderName = '/Volumes/dataFast/toConvert/';
    experimentName = '2018_11_01_looming';%'2018_07_19_hungry';
    datafolder = [baseFolderName,experimentName,'/tiff_stacks/'];     % Folder containing folder of tiff stacks
    savePath = [datafolder,'_matfiles/'];
    mkdir(savePath);
else
    addpath(genpath('C:\Users\Axel-SCAPE\Documents\_code\analysis\NoRMCorre'))
    baseFolderName = 'D:\';%'E:\SCAPE_DATA\';%'W:\scratchData\'; %'D:\';
    experimentName = '2019_04_22_Nsyb_NLS6s_Su_walk';%'2018_11_29_looming';%'2018_07_19_hungry';
    datafolder = [baseFolderName,experimentName,'\tiff_stacks\'];%[baseFolderName,experimentName,'\tiff_stacks\'];     % Folder containing folder of tiff stacks
    savePath = [datafolder,'_matfiles_new\'];
    mkdir(savePath);
end

chunkSize = 201;%201; %871;%677;%751;                                                    % number of time steps per matfile % make it bigger (around 300)
zeroThreshold = 0; %120;                                                % threshold for re-zeroing data % don't change
timeToExclude = 420; %360;%180; %360; %450; %360; %480;%240; %390; %300;%180;%210; %180;%150;%180;%1; %180;%300; %                                               % number of time steps at beginning to skip %framerate * 30 seconds
cropping.x = 13:106; %32:115;%35:126;%40:188;%35:127;%37:124;%50:119; %39:104;%34:95;%35:126;%17:91; %22:143; %38:118;%32:89;%41:120;%28:104;%30:115;%27:129; %38:141; %46:154;%37:130; %28:132;%32:135;%18:138;%43:174; %42:172;%49:142; %58:175; %57:170;%50:170;%68:177;%32:164;%25:180;%40:195; %25:145; %43:115;%36:119;%39:167; 10/31/2018 fly1 IRtest
cropping.y = 54:332;%41:297;%28:329; %51:346;%53:333;%54:324;%80:322;%78:348;%90:335;%41:322;%78:332;%71:334; %75:354;%70:308;%67:354;%21:310;%23:330;%71:324; %59:334; %51:329;%29:355; %20:275;%39:322; %23:320;%62:340; %16:283;%37:257;%18:280; %35:315;%30:330;%53:343;%58:339;%25:360;%75:350; %60:333; %590:862;%12:290;%[];
cropping.z = 1:156; %11:177;%58:215; %3:177;%9:193;%37:214;%28:150; %38:195;%28:129; %45:196; %7:176; %8:181; %44:179;%38:102;%19:179;%12:154;%60:156;%19:192; %17:214; %19:150; %16:174;%7:207; %20:195;%4:214;%33:256;%47:218;%32:296; %45:180; %56:255;%40:260;%1:295;%42:246;%1:204;%12:190; %7:180; %13:149;%42:184;%29:149;

d = dir([datafolder,'fly1*']);                                        % list of datasets (folders)
templateID = ones(size(d));                                           % specify where to get the template. set as 1:length(d) to give each experiment its own template.  Set as ones(size(d)) to use a single template

%m = matfile('E:\SCAPE_DATA\2018_10_25_IRtest\tiff_stacks\_matfiles_new\2018_10_25_IRtest_G_fly1_run1_181.mat');
%template = m.template;
template = [];

for j=1:length(d)
    disp(j)
    trialfolder = [datafolder,d(j).name,'\'];
    if j == 1
        makeManyMatfilesFromTiffStack_dualColor(trialfolder,savePath,experimentName,chunkSize,zeroThreshold,timeToExclude,cropping,template);
    else
        makeManyMatfilesFromTiffStack_dualColor(trialfolder,savePath,experimentName,chunkSize,zeroThreshold,timeToExclude,cropping,template);
    end
end