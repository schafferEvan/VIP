
%experimentFolder = '/Volumes/SCAPEdata1/scratchData/2018_11_01_looming/tiff_stacks/_matfiles_pt1/'; %'/Volumes/data/toRegister/_matfiles/';
experimentFolder = 'D:\2019_04_22_Nsyb_NLS6s_Su_walk\tiff_stacks\_matfiles_new\'; %'D:\2019_03_27_Nsyb_NLS6s_Su_feast\tiff_stacks\_matfiles_new\fly3\';%'/Volumes/SCAPEdata1/scratchData/2018_12_11_hunger/_matfiles_new/'; %'/Volumes/SCAPEdata1/scratchData/2018_11_28_looming/_matfiles_new/'; %'/Volumes/dataFast/toRegister/tmp/'; %2018_11_01_looming/'; %'/Volumes/data/toRegister/_matfiles/';

trialsTmp = dir([experimentFolder,'*.mat']);
isRealFile = zeros(size(trialsTmp));
for j=1:length(trialsTmp)
    isRealFile(j) = ~strcmp(trialsTmp(j).name(1),'.');
end

idx = find(isRealFile);
trials = trialsTmp(idx(1));
for j=2:length(idx)
    trials(j) = trialsTmp(idx(j));
end

frameNum = zeros(size(trials));
runNum = zeros(size(trials));
for j=1:length(trials)
    runLoc = strfind(trials(j).name,'run');
    underscoreLoc = strfind(trials(j).name,'_');
    if ~isempty(runLoc)
        uRel = find(underscoreLoc>runLoc,1,'first');
        runNum(j) = str2double(trials(j).name(runLoc+3:underscoreLoc(uRel)-1));
    else
        runNum(j) = 1;
    end
    frameNum(j) = str2double(trials(j).name(underscoreLoc(end)+1:(end-4)));
    
end
disp(frameNum)
disp(runNum)
[~,trialOrder] = sort( runNum*10^ceil(log10(max(frameNum))) + frameNum, 'ascend' );

%xRange = 30:130;
load([experimentFolder,trials(1).name],'Y');
sz = size(Y);
stepSize = length(trials);
md = round(length(trials)/2); %ceil(length(trials)/2);
m = matfile([experimentFolder,trials(trialOrder(md)).name],'Writable',true);
commonTemp = m.Y(:,:,:,1);
templates = repmat(commonTemp,1,1,1,2); %zeros(length(xRange),sz(2),sz(3),floor(md/stepSize));

for i=md:-1:1
    disp(i)
    m = matfile([experimentFolder,trials(trialOrder(i)).name],'Writable',true);
    y = m.Y(:,:,:,1);
    templates(:,:,:,1) = y; % the first template is always local
    m.templates = templates;
end
for i=md+1:length(trials)
    disp(i)
    m = matfile([experimentFolder,trials(trialOrder(i)).name],'Writable',true);
    y = m.Y(:,:,:,1);
    templates(:,:,:,1) = y; % the first template is always local
    m.templates = templates;
end