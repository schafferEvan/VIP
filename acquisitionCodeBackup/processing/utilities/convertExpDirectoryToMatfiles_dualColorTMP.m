
% reads a directory of folders, each containing tiff stacks, and exports a
% set of matfiles containing chunks of data and a template


%pause(8*60*60);


if ismac
    addpath(genpath('/Users/evan/Dropbox/_code/NoRMCorre'))
    baseFolderName = '/Volumes/dataFast/toConvert/';
    experimentName = '2018_11_01_looming';%'2018_07_19_hungry';
    datafolder = [baseFolderName,experimentName,'/tiff_stacks/'];     % Folder containing folder of tiff stacks
    savePath = [datafolder,'_matfiles/'];
    mkdir(savePath);
else
    addpath(genpath('C:\Users\Axel-SCAPE\Documents\_code\analysis\NoRMCorre'))
    baseFolderName = 'E:\SCAPE_DATA\';%'W:\scratchData\'; %'E:\SCAPE_DATA\';%'D:\';
    experimentName = '2018_11_01_looming';%'2018_10_31_IRtest';%'2018_07_19_hungry';
    datafolder = [baseFolderName,experimentName,'\tiff_stacks\'];     % Folder containing folder of tiff stacks
    savePath = [datafolder,'_matfiles_new\'];
    mkdir(savePath);
end

chunkSize = 201;%201; %871;%677;%751;                                                    % number of time steps per matfile % make it bigger (around 300)
zeroThreshold = 0; %120;                                                % threshold for re-zeroing data % don't change
timeToExclude = 150;%150;%180;%1; %180;%300; %                                               % number of time steps at beginning to skip %framerate * 30 seconds
cropping.x = 57:170;%50:170;%68:177;%32:164;%25:180;%40:195; %25:145; %43:115;%36:119;%39:167; 10/31/2018 fly1 IRtest
cropping.y = 35:315;%30:330;%53:343;%58:339;%25:360;%75:350; %60:333; %590:862;%12:290;%[];
cropping.z = 56:255;%40:260;%1:295;%42:246;%1:204;%12:190; %7:180; %13:149;%42:184;%29:149;

d = dir([datafolder,'fly1_run*']);                                        % list of datasets (folders)
templateID = ones(size(d));                                           % specify where to get the template. set as 1:length(d) to give each experiment its own template.  Set as ones(size(d)) to use a single template

%m = matfile('E:\SCAPE_DATA\2018_10_25_IRtest\tiff_stacks\_matfiles_new\2018_10_25_IRtest_G_fly1_run1_181.mat');
%template = m.template;

for j=1:length(d)
    trialfolder = [datafolder,d(j).name,'\'];
     if j == 1
         template = makeManyMatfilesFromTiffStack_dualColor(trialfolder,savePath,experimentName,chunkSize,zeroThreshold,timeToExclude,cropping);
     else
        makeManyMatfilesFromTiffStack_dualColor(trialfolder,savePath,experimentName,chunkSize,zeroThreshold,timeToExclude,cropping,template);
     end
end