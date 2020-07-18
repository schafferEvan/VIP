function [im_bw_out,cc,regionProps] = cell_segment_watershed(im, confocal_40x)
% this function takes as an input an image with bright blobs and segments
% it using a watershed filter.

addpath(genpath('..'))

% Set parameters
if nargin<2
    confocal_40x = false;
end
if confocal_40x
    % parameters for high-res confocal
    landmarkParams.s1 = 1.5; %3.0;    %0.8;    % sigma of smaller Gaussian
    landmarkParams.s2 = 10; %16;     %2.0;    % sigma of larger Gaussian
    minObjSize=100; %750;             % 27;              % bound on min bump size to pass mask
    minObjDims=[6,6,4];       % [4,4,4];         % bound on min extent of a cell in x,y,z
else
    % standard parameters for scape experiment
    landmarkParams.s1 = 0.8;    % sigma of smaller Gaussian
    landmarkParams.s2 = 2.0;    % sigma of larger Gaussian
    minObjSize= 27;              % bound on min bump size to pass mask
    minObjDims= [4,4,4];         % bound on min extent of a cell in x,y,z
end

% DoG filter to enhance contrast
[~,im_DoG] = findLandmarkPointDoG(im,landmarkParams);
im_DoG_norm=normalizeRangeES(im_DoG);

% define threshold for signal vs noise using GMM
% M=max(im_DoG_norm,[],3);M=M(:);                 % dorsal MIP
M = im_DoG_norm(:);
s = [0,0];
GMModel = fitgmdist(M,2);                       % Gaussian mixture model of of MIP image
P = posterior(GMModel,M);                       % Posterior
s(1) = std( M( P(:,1) > P(:,2) ) );             % std of voxels in dist 1
s(2) = std( M( P(:,1) < P(:,2) ) );             % std of voxels in dist 2
if s(1)<s(2);   noise_id = [true,false];
else;           noise_id = [false,true];
end
is_brain_thresh = GMModel.mu(noise_id) + 5*s(noise_id);   % thresh is 5 sigma above noise mean


% create mask to evaluate Hessian where F is large, ignoring small maxima
mask_bw=im_DoG_norm>is_brain_thresh;
mask_bw=AreaFilter(mask_bw,minObjSize,[],6);


% watershed filter to define cell bounds
watershed_mask=watershed(-im_DoG,6);            % choose connectivity of basins (6,18,26) 

% basins of attraction, labeled with individual identity for visibility
watershed_im = (im_DoG>0).*double(watershed_mask).*mask_bw;
im_bw = watershed_im>0;

%remove small objects by volume
im_bw_out=AreaFilter(im_bw,minObjSize,[],6);

%remove small objects by extent in x,y,z
im_bw_out=extentFilter(im_bw_out,minObjDims);

%for confocal (high res), remove non-spherical objects
if confocal_40x
    im_bw_out=sphericityFilter(im_bw_out);
    im_bw_out=maxFilter(im_bw_out, im_DoG);
end

%compile results of all sub images into final segmented mask.
regionProps = getRegionProps(im_bw_out);
% figure; imagesc(squeeze(sum(im_bw_out,3)))
disp(['found ',num2str(length(regionProps.vol)),' objects'])
cc=regionProps.cc;




function Imout=sphericityFilter(Imin)
% generate image without ROIs that fail any minObjDims
sphericity_th = 0.4;
regionProps = getRegionProps(Imin);
cc=bwconncomp(Imin,6);
sok = zeros(length(regionProps.blobStats),1);
for j=1:length(regionProps.blobStats)
    sok(j) = regionProps.sphericity(j)>sphericity_th;
end
cc.PixelIdxList(~sok)=[];
Imout=zeros(size(Imin));
Imout(cell2mat(cc.PixelIdxList'))=1;


function Imout=maxFilter(Imin, analog_im)
% generate image without ROIs that fail any minObjDims
max_th = 0.05*max(analog_im(:));
regionProps = getRegionProps(Imin);
cc=bwconncomp(Imin,6);
sok = zeros(length(regionProps.blobStats),1);
for j=1:length(regionProps.blobStats)
    sok(j) = mean(analog_im(cc.PixelIdxList{j}))>max_th;
end
cc.PixelIdxList(~sok)=[];
Imout=zeros(size(Imin));
Imout(cell2mat(cc.PixelIdxList'))=1;


function Imout=extentFilter(Imin,minObjDims)
% generate image without ROIs that fail any minObjDims
regionProps = getRegionProps(Imin);
cc=bwconncomp(Imin,6);

xok = zeros(length(regionProps.blobStats),1);
yok = zeros(length(regionProps.blobStats),1);
zok = zeros(length(regionProps.blobStats),1);
for j=1:length(regionProps.blobStats)
    xok(j) = regionProps.blobStats(j).BoundingBox(4)>=minObjDims(1);
    yok(j) = regionProps.blobStats(j).BoundingBox(5)>=minObjDims(2);
    zok(j) = regionProps.blobStats(j).BoundingBox(6)>=minObjDims(3);
end
cc.PixelIdxList(~xok | ~yok | ~zok)=[];
Imout=zeros(size(Imin));
Imout(cell2mat(cc.PixelIdxList'))=1;



function regionProps = getRegionProps(BW)
cc=bwconncomp(BW,6);
blobStats=regionprops(cc,'Area','BoundingBox','Centroid');

% calculate blob parameters
perim_im=bwperim(BW,26);
vol=cellfun(@(x) length(x),cc.PixelIdxList);
perim=cellfun(@(x) sum(perim_im(x)),cc.PixelIdxList);
sphericity=((pi*36*vol.^2).^(1/3)./perim);

regionProps.cc = cc;
regionProps.vol = vol;
regionProps.sphericity = sphericity;
regionProps.blobStats = blobStats;
figure; plot(vol,sphericity,'o')


function V=normalizeRangeES(X)
%The function takes a scalar or nd matrix X and normalize the matrix to a
%minimum of 0 and a maximum of 1. If all the values in the matrix are the
%same, they are all replaced with ones.
if ~isempty(X)
    if max(X(:))==min(X(:))
        V=ones(size(X));
    else
        V=(X-min(X(:)))/(max(X(:))-min(X(:)));
    end
else
    V=X;
end

