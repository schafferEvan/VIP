
experimentFolder = 'D:\2019_03_26_Nsyb_NLS6s_Su_feast\tiff_stacks\_matfiles_new\';

trials = dir([experimentFolder,'*.mat']);
frameNum = zeros(size(trials));
runNum = zeros(size(trials));
for j=1:length(trials)
    runLoc = strfind(trials(j).name,'run');
    underscoreLoc = strfind(trials(j).name,'_');
    frameNum(j) = str2double(trials(j).name(underscoreLoc(end)+1:(end-4)));
    runNum(j) = str2double(trials(j).name(runLoc+3:underscoreLoc(end-1)-1));
end
[~,trialOrder] = sort(runNum*10^ceil(log10(max(frameNum))) + frameNum, 'ascend' );

%xRange = 30:130;
load([experimentFolder,trials(1).name],'Y');
sz = size(Y);
stepSize = 3;
md = round(length(trials)/2);
m = matfile([experimentFolder,trials(trialOrder(md)).name],'Writable',true);
commonTemp = m.Y(:,:,:,1);
SharedTemplates = repmat(commonTemp,1,1,1,length(1:ceil(md/stepSize)+1)); %zeros(length(xRange),sz(2),sz(3),floor(md/stepSize));


ctr=0;
for i=md:-1:1
    disp(i)
    m = matfile([experimentFolder,trials(trialOrder(i)).name],'Writable',true);
    y = m.Y(:,:,:,1);
    if ~mod(md-i,stepSize)
        SharedTemplates(:,:,:,1:floor(md/stepSize)-ctr) = repmat(y,1,1,1,length(1:floor(md/stepSize)-ctr));
        ctr = ctr+1;
        templates = SharedTemplates(:,:,:,end-ctr:end);
    else
        templates = SharedTemplates(:,:,:,end-ctr:end);
        templates = cat(4,y,templates); % the first template is always local
    end
    m.templates = templates;
end

ctr=0;
for i=md+1:length(trials)
    disp(i)
    m = matfile([experimentFolder,trials(trialOrder(i)).name],'Writable',true);
    y = m.Y(:,:,:,1);
    if ~mod(md-i,stepSize)
        SharedTemplates(:,:,:,1:floor(md/stepSize)-ctr) = repmat(y,1,1,1,length(1:floor(md/stepSize)-ctr));
        ctr = ctr+1;
        templates = SharedTemplates(:,:,:,end-ctr:end);
    else
        templates = SharedTemplates(:,:,:,end-ctr:end);
        templates = cat(4,y,templates); % the first template is always local
    end
    m.templates = templates;
end
