
function GMMreg_toCommonCoords(codePath, experimentFolder, figureFolder)

addpath(genpath(codePath))

% raw pointset to be aligned
rawFile = load([experimentFolder,'cc.mat']);
d = dir([experimentFolder,'*Agood.mat']);
goodIds_raw = load([experimentFolder,d.name],'goodIds');
goodIds_raw = find(goodIds_raw.goodIds);
raw = zeros(length(goodIds_raw),3);
for k=1:size(raw,1)
    raw(k,:) = rawFile.regionProps.blobStats(goodIds_raw(k)).Centroid;
end

% hardcoded refFile is constant for all datasets
disp('-- aligning to dorsoposterior template --')
load('/Volumes/SCAPEdata1/finalData/_templates/dorsopost_25deg_template.mat','ref');

% disp('-- aligning to 2018/08/24 NsybNLS fly2 --')
% refFile = load('/Volumes/SCAPEdata1/finalData/2018_08_24_NsybNLS_odors/fly2/Yproj/cc.mat');
% goodIds_ref = load('/Volumes/SCAPEdata1/finalData/2018_08_24_NsybNLS_odors/fly2/Yproj/2018_08_24_NsybNLS_odors_fly2_Agood.mat','goodIds');
% goodIds_ref = find(goodIds_ref.goodIds);
% ref = zeros(length(goodIds_ref),3);
% for k=1:size(ref,1)
%     ref(k,:) = refFile.regionProps.blobStats(goodIds_ref(k)).Centroid;
% end
        
% parameters
motion = 'tps';         % motion model ('tps','affine3d')
scale = 50;             % variance around each point
interval = 5; %10;      % size of TPS grid

if contains(experimentFolder, '2018_08_24_NsybNLS_odors/fly2')
    % skip alignment for template dataset
    aligned = ref;
else
    % align point sets using TPS GMMreg
    [config] = initialize_config_extended(raw, ref, motion, scale, raw, ref, interval);
    [~, aligned] =  gmmreg_L2(config);
end

save([experimentFolder,'registered_pointset.mat'],'aligned','raw')

% visualize results
raw(:,3)=-raw(:,3);
ref(:,3)=-ref(:,3);
aligned(:,3)=-aligned(:,3);
DisplayPoints3DPretty(raw, ref, aligned);

% make movie
if ~exist(figureFolder,'dir'); mkdir(figureFolder); end
OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
CaptureFigVid_multiAxis([-20,40;-110,40;-190,89;-290,40;-380,40], [figureFolder,'registered_pointset.mp4'],OptionZ)



