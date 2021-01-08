% Scripts for running post-processing on Terremoto

codePath = '../../';
experimentFolder = '/moto/axs/projects/registered/2018_08_24_fly2/fly2/';

% function run_watershed_segmentation(codePath, experimentFolder)

addpath(genpath(codePath))
    
load([experimentFolder,'Yproj/Ysum.mat'],'Ysum','Rsum')

[~,Ycc,regionProps_green] = cell_segment_watershed(Ysum);
% [im_bw_out,Ycc,regionProps] = segmentHessian_SCAPE_new(Ysum);
if ~isempty(Rsum)
    [~,Rcc,regionProps] = cell_segment_watershed(Rsum);
else
    Rcc = [];
end
save([experimentFolder,'Yproj/cc.mat'],'Ycc','Rcc','regionProps','regionProps_green')