
function make_quickMovie_twoColor(codePath, experimentFolder, savePath)
% skips frames to quickly make movie containing sample frames from entire
% experiment (primarily to check motion correction)

addpath(genpath(codePath))
infoFile = dir([experimentFolder,'info/*.mat']);

trials = sortExperimentDirectory(experimentFolder,'reg');
params.trialName = [trials(end).name(1:end-4),'_small.avi']; %[trials(end).name(1:runLoc-1),'all.avi'];
params.savePath = savePath;
if ~isfolder(savePath); mkdir(savePath); end

params.concatenate = true;
j=1;
while j<=length(infoFile)
    try
        load([experimentFolder,'info/',infoFile(j).name]);
    catch
        j=j+1;
    end
end
params.Ttot = 0;
params.acqRate = round(info.daq.scanRate); %10; % volumes per second

stepSize = 50;

% select maxY to set color scale
m = matfile([experimentFolder,trials(2).name]);
sz = size(m,'Y');
Y = single( m.Y(:,:,:,sz(4)-50:sz(4) ) );
params.maxY = max(Y(:));

try
    params.xum = sz(1)*info.GUIcalFactors.x_umPerPix;
catch
    params.xum = sz(1)*info.GUIcalFactors.xK_umPerVolt*info.daq.scanAngle/(info.daq.pixelsPerLine-1);
end
params.yum = sz(2)*info.GUIcalFactors.y_umPerPix;
params.zum = sz(3)*info.GUIcalFactors.z_umPerPix;

rTh = 250;
mask = m.R(:,:,:,1)>rTh;


disp(['making quick movie of ',num2str(length(trials)),' trials'])
params.kTot=0;
%% load in registered data and find ROIs
for i=1:length(trials)
    params.Ttot = 0;
    
    trialPath = [experimentFolder,trials(i).name];
    display(['file ',num2str(i),': ',trials(i).name])

    m = matfile(trialPath);
    try
        sz = size(m,'Ysmall');
    catch
        sz = size(m,'Y');
    end
    pts = 1:stepSize:sz(end);
    if length(pts)<2; pts=[1,2]; end
    Y = single( m.Y(:,:,:,pts) );
    szy = size(Y);
    Y = Y/params.maxY;
    
    stimOn = [inf,inf];
    params.vid = plot4Dproj_ES(Y, mask, [szy(1),szy(2),szy(3)],params,stimOn);
    params.Ttot = params.Ttot + .1*round(10*szy(4)/params.acqRate);
    params.kTot=params.kTot+szy(4);
    
end
close(params.vid)