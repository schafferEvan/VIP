
function GMMreg_toCommonCoords(codePath, experimentFolder, figureFolder, align_on_green)

addpath(genpath(codePath))
% rawFile = load([experimentFolder,'cc.mat']);
rawFile = load([experimentFolder,'centroids.mat']);

if nargin<4
    align_on_green=0;
else
    align_on_green = eval(align_on_green);
end

% raw pointset to be aligned
if align_on_green
    % align green centroids (for sparse red pan green datasets)
    %     raw = zeros(length(rawFile.regionProps_green.blobStats),3);
    %     for k=1:size(raw,1)
    %         raw(k,:) = rawFile.regionProps_green.blobStats(k).Centroid;
    %     end
    raw = rawFile.centroids_fromGreen;
    
    % append red centroids
    d = dir([experimentFolder,'*Agood.mat']);
    goodIds_raw = load([experimentFolder,d.name],'goodIds');
    %     goodIds_raw = find(goodIds_raw.goodIds);
    %     raw_red = zeros(length(goodIds_raw),3);
    %     for k=1:size(raw_red,1)
    %         raw_red(k,:) = rawFile.regionProps.blobStats(goodIds_raw(k)).Centroid;
    %     end
    raw_red = rawFile.centroids(goodIds_raw.goodIds,:);
    is_red = logical([zeros(size(raw,1),1);ones(size(raw_red,1),1)]);
    raw = [raw;raw_red];
else
    % default: align red centroids
    d = dir([experimentFolder,'*Agood.mat']);
    goodIds_raw = load([experimentFolder,d.name],'goodIds');
    %     goodIds_raw = find(goodIds_raw.goodIds);
    %     raw = zeros(length(goodIds_raw),3);
    %     for k=1:size(raw,1)
    %         raw(k,:) = rawFile.regionProps.blobStats(goodIds_raw(k)).Centroid;
    %     end
    is_red = true(sum(goodIds_raw.goodIds),1);
    raw = rawFile.centroids(goodIds_raw.goodIds,:);
end

% hardcoded refFile is constant for all datasets
disp('-- aligning to dorsoposterior template --')
try
    load('/Volumes/SCAPEdata1/finalData/_templates/dorsopost_template.mat','ref');
catch
    disp('referencing experiment folder')
    s=strfind(experimentFolder,'/');
    templatePath=experimentFolder(1:s(end-3));
    load([templatePath,'_templates/dorsopost_template.mat'],'ref');
end

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
interval = 5; %10; %5;      % size of TPS grid

% if contains(experimentFolder, '2018_08_24_NsybNLS_odors/fly2')
%     % skip alignment for template dataset
%     aligned = ref;
% else

% align point sets using TPS GMMreg
[config] = initialize_config_extended(raw, ref, motion, scale, raw, ref, interval);
[~, aligned] =  gmmreg_L2(config);
% end

if align_on_green
    save([experimentFolder,'registered_pointset.mat'],'aligned','raw','is_red')
else
    save([experimentFolder,'registered_pointset.mat'],'aligned','raw')
end
% visualize results
raw(:,3)=-raw(:,3);
ref(:,3)=-ref(:,3);
aligned(:,3)=-aligned(:,3);
DisplayPoints3DPretty(raw, ref, aligned, is_red);

% make movie
if ~exist(figureFolder,'dir'); mkdir(figureFolder); end
OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
CaptureFigVid_multiAxis([-20,40;-110,40;-190,89;-290,40;-380,40], [figureFolder,'registered_pointset.mp4'],OptionZ)
disp('--- COMPLETED ALIGNMENT OF POINT SET TO REFERENCE IMAGE ---')


