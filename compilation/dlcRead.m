
function dlcArray = dlcRead(filename,npts)
% Import data from deeplabcut csv file.

if ~nargin
    filename = '/Volumes/SCAPEdata1/finalData/2019_07_01_Nsyb_NLS6s_walk/fly2/Yproj/701_fly2DeepCut_resnet50_RunningJul3shuffle1_1030000.csv';
    npts=8; %number of tracked points
    %filename = '/Volumes/SCAPEdata1/finalData/2019_07_01_Nsyb_NLS6s_walk/fly2/Yproj/701_fly2DeepCut_resnet50_RunningJul3shuffle1_1030000.csv';
end
delimiter = ',';
startRow = 4;

% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

% Close the text file.
fclose(fileID);

% Create output variable
l = {'x','y','likelihood'};
labels = genvarname( repmat(l,1,npts), l );
dlcArray = table(dataArray{2:end-1}, 'VariableNames', labels);
