function [trials, trialOrder, runNum, runNumConvert, frameNum] = sortExperimentDirectory(experimentFolder)
% File sorting in experimental directory.

trials = dir([experimentFolder,'*.mat']); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%dir([experimentFolder,'*.h5']);
frameNum = zeros(size(trials));
runNum = cell(size(trials));
runNumConvert = zeros(size(trials));
for j=1:length(trials)
    m = matfile([experimentFolder,trials(j).name]); % this is to check if any file is corrupted
    runLoc = strfind(trials(j).name,'run');
    underscoreLoc = strfind(trials(j).name,'_');
    regLoc = strfind(trials(j).name,'reg');
    if isempty(regLoc)
        regLoc = strfind(trials(j).name,'.');
    end
    if ~isempty(runLoc)
        runNum{j} = trials(j).name(runLoc+3:underscoreLoc(end)-1);
    else
        runNum{j} = '1';
    end
    frameNum(j) = str2double(trials(j).name(underscoreLoc(end)+1:regLoc(1)-1));
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

[~,trialOrder] = sort( runNumConvert*10^(1+ceil(log10(max(frameNum)))) + frameNum, 'ascend' );
end

