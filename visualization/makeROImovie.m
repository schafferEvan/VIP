
function makeROImovie(raw_image, footprints, save_path, info, is_confocal)
% addpath(genpath( '~/Dropbox/_code/'))
% expDir = '/Users/evan/Dropbox/_sandbox/sourceExtraction/good/';
% expID = '0824_f3r1'; %'0312_f4'; %'0327_f3'; %'0824_f2r2_cur'; %'0221'; %'0327_f3'; %'0321'; %
% 
% fromGreenCC = false;
% if fromGreenCC
%     load([expDir,expID,'/post_fromYcc.mat'])
%     load([expDir,expID,'/F.mat'])
%     matfileRaw = matfile([expDir,expID,'/F.mat']);
%     matfilePost = matfile([expDir,expID,'/post_fromYcc.mat']);
%     ccFlag = ' from Ycc';
% else
%     load([expDir,expID,'/post_fromRcc.mat'])
%     load([expDir,expID,'/F_fromRed.mat'])
%     matfileRaw = matfile([expDir,expID,'/F_fromRed.mat']);
%     matfilePost = matfile([expDir,expID,'/post_fromRcc.mat']);
%     ccFlag = ' from Rcc';
% end
% load([expDir,expID,'/Ysum.mat'])
% matfileBeh = matfile([expDir,expID,'/alignedBehavAndStim.mat']);
% infoFile = dir([expDir,expID,'/info/*.mat']);%'/Volumes/dataFast/habaRegistered/2018_08_24_odor/mats/new/fly3run2/'; %'/Users/evan/Desktop/hungerRaw/';%'/Volumes/SCAPEdata1/scratchData/2018_08_01_IRtest/matfiles/registered/';%'/Volumes/data/_scape/data/_outMats/'; %
% if length(infoFile)>1; infoFile=infoFile(end); end
% load([expDir,expID,'/info/',infoFile.name]);
sz = size(raw_image);

A = sparse(prod(footprints.ImageSize), footprints.NumObjects);
mx = zeros(footprints.NumObjects,1);
for j=1:footprints.NumObjects
    A(footprints.PixelIdxList{j},j)=1;
    mx(j) = max(raw_image(footprints.PixelIdxList{j}));
end
num_comps = footprints.NumObjects;




% ------------ DIMENSIONS -----------------------------
if is_confocal
    params.xum = sz(1);
    params.yum = sz(2);
    params.zum = sz(3);
    totWidth = 1.1*params.yum; %(params.zum + params.yum);
    params.zpx = params.zum/totWidth;
    params.ypx = params.yum/totWidth;
    params.xpx = params.xum/totWidth;
else
    try
        params.xum = sz(1)*info.GUIcalFactors.x_umPerPix;
    catch
        %warning('no x_um field available')
        params.xum = sz(1)*info.GUIcalFactors.xK_umPerVolt*info.daq.scanAngle/(info.daq.pixelsPerLine-1);
    end
    params.yum = sz(2)*info.GUIcalFactors.y_umPerPix;
    params.zum = sz(3)*info.GUIcalFactors.z_umPerPix;
    totWidth = 1.1*params.yum; %(params.zum + params.yum);
    params.zpx = params.zum/totWidth;
    params.ypx = params.yum/totWidth;
    params.xpx = params.xum/totWidth;
end


% ------------ MAKE FINDER PLOT OR FINAL PLOT -----------------------------
bkTh = 0.1;
colors = hsv(size(A,2));
colors = colors( randperm(size(A,2)), : );

Ar1 = reshape( full(sum(A*sparse(diag(colors(:,1)')),2)), sz(1),sz(2),sz(3));
Ar2 = reshape( full(sum(A*sparse(diag(colors(:,2)')),2)), sz(1),sz(2),sz(3));
Ar3 = reshape( full(sum(A*sparse(diag(colors(:,3)')),2)), sz(1),sz(2),sz(3));
mx = max([max(Ar1(:)),max(Ar2(:)),max(Ar3(:))]);

vid = VideoWriter([save_path,'_ROIvid.avi']);
vid.FrameRate = 30;
open(vid)

pi.totScale = 0.7; % ratio of height to width for the whole figure. image dims are compensated to keep correct scaling
figw = 800;
pi.figpos = [368 286 figw figw*pi.totScale]; % figure position
f2 = figure('position', pi.figpos);
set(f2,'color','k')
axes('position',[.05  .05/pi.totScale  params.ypx  params.xpx/pi.totScale]); %subplot(1, 3, 1);

disp(['Making movie of ',num2str(size(A,2)),' components'])
for z=1:sz(3)-3
    R = Ar1(:,:,z);
    G = Ar2(:,:,z);
    B = Ar3(:,:,z);
    
    R = R/mx; G = G/mx; B = B/mx;
    RGB = cat(3,R,G,B);
    bk1 = max(R+G+B,[],3); bk1 = bk1/quantile(bk1(:),0.95); bk1(bk1>1)=1;
    al = bk1>bkTh; %squeeze(bk1.*bk.*(bk2>bkTh)); %
    b = squeeze(raw_image(:,:,z));
    b = b/max(raw_image(:));
    b = cat(3,b,b,b);
    Kadj = imadjust(b, [0 1], [0 1], 0.6 );
    
    annotation('line',.875+[0 25*params.ypx/params.yum],[.11 .11],'linewidth',2,'color',[.6 .6 .6])
    annotation('textbox',[.875 .1 25*params.ypx/params.yum .01],'String','25um','color',[.5 .5 .5],'fontsize',12,'linestyle','none','HorizontalAlignment','center')
    annotation('textbox',[.8 .9 .2 .05],'String',[num2str(num_comps),' comps'],'color',[.5 .5 .5],'fontsize',16,'linestyle','none','HorizontalAlignment','center')
    
    imagesc(Kadj)
    rgb_adj = imadjust(RGB, [0.1 0.7], [0 1], 0.6 );
    hold on; im=imagesc(rgb_adj);%/q);
    set(im, 'AlphaData', 4*al);%4*al);
    axis off;
    drawnow
    writeVideo(vid, getframe(f2))

end
close(vid)
