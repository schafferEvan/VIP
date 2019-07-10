%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Volumes/SCAPEdata1/finalData/2019_07_01_Nsyb_NLS6s_walk/fly2/Yproj/701_fly2DeepCut_resnet50_RunningJul3shuffle1_1030000.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2019/07/10 18:38:53

%% Initialize variables.
filename = '/Volumes/SCAPEdata1/finalData/2019_07_01_Nsyb_NLS6s_walk/fly2/Yproj/701_fly2DeepCut_resnet50_RunningJul3shuffle1_1030000.csv';
delimiter = ',';
startRow = 4;

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
fly2DeepCutresnet50RunningJul3shuffle11030000 = table(dataArray{1:end-1}, 'VariableNames', {'scorer','DeepCut_resnet50_RunningJul3shuffle1_1030000','DeepCut_resnet50_RunningJul3shuffle1_1','DeepCut_resnet50_RunningJul3shuffle1_2','DeepCut_resnet50_RunningJul3shuffle1_3','DeepCut_resnet50_RunningJul3shuffle1_4','DeepCut_resnet50_RunningJul3shuffle1_5','DeepCut_resnet50_RunningJul3shuffle1_6','DeepCut_resnet50_RunningJul3shuffle1_7','DeepCut_resnet50_RunningJul3shuffle1_8','DeepCut_resnet50_RunningJul3shuffle1_9','DeepCut_resnet50_RunningJul3shuffle1_10','DeepCut_resnet50_RunningJul3shuffle1_11','DeepCut_resnet50_RunningJul3shuffle1_12','DeepCut_resnet50_RunningJul3shuffle1_13','DeepCut_resnet50_RunningJul3shuffle1_14','DeepCut_resnet50_RunningJul3shuffle1_15','DeepCut_resnet50_RunningJul3shuffle1_16','DeepCut_resnet50_RunningJul3shuffle1_17','DeepCut_resnet50_RunningJul3shuffle1_18','DeepCut_resnet50_RunningJul3shuffle1_19','DeepCut_resnet50_RunningJul3shuffle1_20','DeepCut_resnet50_RunningJul3shuffle1_21','DeepCut_resnet50_RunningJul3shuffle1_22','DeepCut_resnet50_RunningJul3shuffle1_23'});

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;