function checkCorruptedMatfiles(codePath, experimentFolder)
% checks for any corrupted matfiles in directory
 
addpath(genpath(codePath))
trials = sortExperimentDirectory(experimentFolder,'reg');
 
%% load in registered data and find ROIs
for i=1:length(trials)
    trialPath = [experimentFolder,trials(i).name];
 
    try
        m = matfile(trialPath);
    catch
        %warning([trials(i).name, ' is a corrupted file.'])
        disp([trials(i).name, ' cannot be turned into a matfile.'])
    end   
     
    try
        a = load(trialPath);
    catch
        %warning([trials(i).name, ' is a corrupted file.'])
        disp([trials(i).name, ' cannot be loaded.'])
    end  
     
    clear a m
end