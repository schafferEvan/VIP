

function run_watershed_segmentation(codePath, experimentFolder)

addpath(genpath(codePath))
    
load([experimentFolder,'Yproj/Ysum.mat'],'Ysum','Rsum')

[~,Ycc] = cell_segment_watershed(Ysum);
% [im_bw_out,Ycc,regionProps] = segmentHessian_SCAPE_new(Ysum);
if ~isempty(Rsum)
    [~,Rcc,regionProps] = cell_segment_watershed(Rsum);
else
    Rcc = [];
end
save([experimentFolder,'Yproj/cc.mat'],'Ycc','Rcc','regionProps')