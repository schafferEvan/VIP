function makeFootprintFig(fd)
expID=fd.expID;
savePath = fd.savePath;
fnm=fd.fnm;
ccFlag=fd.ccFlag;

pOrd=fd.pOrd;
d1=fd.d1;
d2=fd.d2;
d3=fd.d3;

K=fd.K;
A=fd.A;
goodIds=fd.goodIds;
idx=fd.idx;

Ysum=fd.Ysum;



f2 = figure;
set(gcf,'color','w')
pos=get(gcf,'Position');
pos(3)=800;%pos(3)*1.5;
pos(4)=500;%pos(4)*1.5;
set(gcf,'Position',pos)
f2.InvertHardcopy = 'off';
f2.PaperUnits = 'points';
f2.PaperSize = 1.1*[pos(3) pos(4)];


% make map of footprints --------------------------------------------------
bkTh = 0.1;
groupColors = flipud(jet(K));
R = zeros(d1,d2);
G = zeros(d1,d2);
B = zeros(d1,d2);
Agood = A(:,goodIds);
cIdx = zeros(size(idx));
tmpIdx = 0;

for k=1:K
    kIds = idx==pOrd(k); %goodIds(idx==pOrd(k));
    cIdx(tmpIdx+(1:sum(kIds)))=k; tmpIdx=tmpIdx+sum(kIds);
    Ar = reshape( max(full(Agood(:,kIds)),[],2), d1,d2,d3);
    Ar1 = max(Ar,[],3);
    R = R + Ar1*groupColors(k,1);
    G = G + Ar1*groupColors(k,2);
    B = B + Ar1*groupColors(k,3);
end
%mx = max([max(R(:)),max(G(:)),max(B(:))]);
%R = R/mx; G = G/mx; B = B/mx;
RGB = cat(3,R,G,B);




% show map of footprints --------------------------------------------------
bk1 = max(R+G+B,[],3); bk1 = bk1/quantile(bk1(:),0.95); bk1(bk1>1)=1;
al = bk1>bkTh; %squeeze(bk1.*bk.*(bk2>bkTh)); %
b = squeeze(max(Ysum,[],3));
b = b/max(b(:));
b = cat(3,b,b,b);
Kadj = imadjust(b, [0 1], [0 1], 0.6 );
imagesc(Kadj)
rgb_adj = RGB; %imadjust(RGB, [0.1 0.7], [0 1], 0.6 );
hold on; im=imagesc(rgb_adj);%/q);
set(im, 'AlphaData', 4*al);%4*al);
axis off;
expTitle=strrep(expID,'_','\_');
title([expTitle,ccFlag],'color',[.9,.9,.9]);%,'interpreter','none')

%annotation('line',.905+[0 25*ypx/params.yum],[.05 .05],'linewidth',2,'color',[.6 .6 .6])
%annotation('textbox',[.905 .047 25*ypx/params.yum .01],'String','25um','color',[.5 .5 .5],'fontsize',12,'linestyle','none','HorizontalAlignment','center')
setFigColors;
f2.InvertHardcopy = 'off';
if ~isfolder([savePath,'_plots/']); mkdir([savePath,'_plots/',fnm]); end
expNameHandle=strrep(expID,'/','_');
saveas(f2, [savePath,'_plots/',expNameHandle,'_ROIs_from',ccFlag(end-2:end),'.png'])



function setFigColors
set(gca,'Fontsize',14)
set(gca,'color','none');%'k')
set(gcf,'color','none');%'k');
set(gca,'xcolor',[.9 .9 .9]);
set(gca,'ycolor',[.9 .9 .9]);