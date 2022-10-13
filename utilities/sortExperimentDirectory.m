function [trials, trialOrder, runNum, runNumConvert, frameNum] = sortExperimentDirectory(experimentFolder, typeFile, sortOutput)
% File sorting in experimental directory.

trials = dir([experimentFolder,'*.mat']); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%dir([experimentFolder,'*.h5']);

pos = [];
for i =1:numel(trials)
    trialName = trials(i).name;
    if trialName(1) == '.'
        pos = [pos, i];
    end
end 

trials(pos) = [];


frameNum = zeros(size(trials));
runNum = cell(size(trials));
runNumConvert = zeros(size(trials));
for j=1:length(trials)
    %m = matfile([experimentFolder,trials(j).name]); % this is to check if any file is corrupted
    runLoc = strfind(trials(j).name,'run');
    underscoreLoc = strfind(trials(j).name,'_');
    typeLoc = strfind(trials(j).name, typeFile);
    %underscoreInType = strfind(typeFile,'_');
    if isempty(typeLoc)
        typeLoc = strfind(trials(j).name,'.');
    end
    if ~isempty(runLoc)
        runNum{j} = trials(j).name(runLoc+3:underscoreLoc(end)-1);
    else
        runNum{j} = '1';
    end
  
    frameNum(j) = str2double(trials(j).name(underscoreLoc(end)+1:typeLoc(1)-1));
     
    uLocTemp = strfind(runNum{j},'_');
    numTemp = runNum{j};
    if length(uLocTemp) > 0
        numTemp(uLocTemp(1)) = '.';
    end    
    if length(uLocTemp) > 1    
        numTemp(uLocTemp(2:end)) = '';
    end
    runNumConvert(j) = str2double(numTemp);
    
end

if strcmp(typeFile, '_info') == 1
    [~,trialOrder] = sort( runNumConvert, 'ascend' );
else    
    [~,trialOrder] = sort( runNumConvert*10^(1+ceil(log10(max(frameNum)))) + frameNum, 'ascend' );
end

if nargin<3
    sortOutput = true;
end

if sortOutput
    trials = trials(trialOrder);
    runNumConvert = runNumConvert(trialOrder); 
end       

