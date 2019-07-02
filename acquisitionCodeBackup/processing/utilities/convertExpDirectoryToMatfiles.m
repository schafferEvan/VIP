
% reads a directory of folders, each containing tiff stacks, and exports a
% set of matfiles containing chunks of data and a template

if ismac
    addpath(genpath('/Users/evan/Dropbox/_code/NoRMCorre'))
    baseFolderName = '/Volumes/dataFast/toConvert/';
    experimentName = '2018_07_26_IRlaser';%'2018_07_19_hungry';
    datafolder = [baseFolderName,experimentName,'/tiff_stacks/'];     % Folder containing folder of tiff stacks
    savePath = [datafolder,'_matfiles/'];
    mkdir(savePath);
else
    baseFolderName = 'D:\'; %'E:\SCAPE_DATA\';%'D:\';
    experimentName ='2019_03_19_ThR_NLS6s_Su';% '2018_10_24_MCFO_IRtest';%'2018_07_19_hungry';
    datafolder = [baseFolderName,experimentName,'\tiff_stacks\'];     % Folder containing folder of tiff stacks
    savePath = [datafolder,'_matfiles\'];
    mkdir(savePath);
end

chunkSize = 300;%677;%751;                                                    % number of time steps per matfile
zeroThreshold = 0;                                                % threshold for re-zeroing data
timeToExclude = 300;%900;%180;%300; %                                               % number of time steps at beginning to skip
cropping.x = 30:127;%225:96;%58:215;%54:175;%1:155;%36:119;%39:167;
cropping.y = 11:378; %81:347;%43:283;%30:310;%9:282;%12:290;%[];
cropping.z = 22:195; %16:169;%54:199;%59:174;%21:187;%42:184;%29:149;

d = dir([datafolder,'fly1*']);                                        % list of datasets (folders)
templateID = ones(size(d));                                           % specify where to get the template. set as 1:length(d) to give each experiment its own template.  Set as ones(size(d)) to use a single template

% m = matfile('E:\SCAPE_DATA\2018_08_16_IRtest\tiff_stacks\_matfiles\2018_08_16_IRtest_fly1_run1_run1_301.mat');
% template = m.template;

for j=1:length(d)
    trialfolder = [datafolder,d(j).name,'\'];
    if j == 1
        template = makeManyMatfilesFromTiffStack(trialfolder,savePath,experimentName,chunkSize,zeroThreshold,timeToExclude,cropping);
    else
        makeManyMatfilesFromTiffStack(trialfolder,savePath,experimentName,chunkSize,zeroThreshold,timeToExclude,cropping,template);
    end
end