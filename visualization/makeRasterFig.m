function makeRasterFig(fd)
expID=fd.expID;
savePath = fd.savePath;
dFlag=fd.dFlag;
fnm=fd.fnm;
ccFlag=fd.ccFlag;

pData=fd.pData;
pOrdFull=fd.pOrdFull;
pOrd=fd.pOrd;
d1=fd.d1;
d2=fd.d2;
d3=fd.d3;
behNorm=fd.behNorm;
time = fd.time;
alignedBehavior=fd.alignedBehavior;
K=fd.K;
A=fd.A;
goodIds=fd.goodIds;
idx=fd.idx;

Ysum=fd.Ysum;



% generate xticks and labels
tks = floor(time(1)/60):floor(time(end)/60);
tkLoc = zeros(size(tks));
tkLab = cell(size(tks));
for j=1:length(tks)
    if mod(tks(j),2); tkLab{j}='';
    else; tkLab{j}=num2str(j);
    end
    [~,tkLoc(j)] = min(abs(time-60*tks(j)));
end

f2 = figure;
set(gcf,'color','w')
pos=get(gcf,'Position');
pos(3)=500;%pos(3)*1.5;
pos(4)=800;%pos(4)*1.5;
set(gcf,'Position',pos)
f2.InvertHardcopy = 'off';
f2.PaperUnits = 'points';
f2.PaperSize = 1.1*[pos(3) pos(4)];

subplot(5,1,2:5); imagesc(pData(pOrdFull,:)); xlim([1 size(pData,2)]); caxis([0.1 0.9])
ylabel('Cell Number','fontsize',12); xlabel('Time (min)','fontsize',12);
set(gca,'XTick',tkLoc); set(gca,'XTickLabel',tkLab);
colormap hot
setFigColors;
tmpA=get(gca,'Position');

subplot(20,1,2); plot(behNorm,'color',[.9,.9,.9],'linewidth',1); xlim([1 size(pData,2)]); ylim([min(behNorm) max(behNorm)]); axis off
tmpB=get(gca,'Position');
set(gca,'Position',[tmpB(1),tmpB(2)-.4*(tmpB(2)-tmpA(2)-tmpA(4)),tmpB(3),tmpB(4)+.4*(tmpB(2)-tmpA(2)-tmpA(4))]);

subplot(20,1,1); plot(alignedBehavior.stim,'linewidth',1); hold all; plot(alignedBehavior.drink,'linewidth',1); xlim([1 size(pData,2)]); ylim([0 1]); axis off
tmpC=get(gca,'Position');
set(gca,'Position',[tmpC(1),tmpC(2)-.4*(tmpC(2)-tmpB(2)-tmpB(4)),tmpC(3),tmpC(4)+.4*(tmpC(2)-tmpB(2)-tmpB(4))]);
expTitle=strrep(expID,'_','\_');
title([expTitle,dFlag,ccFlag],'color',[.9,.9,.9]);%,'interpreter','none')

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



% make colorbar to relate footprints to traces ----------------------------
atmp = tmpA;%get(gca,'position');
axes('position',[atmp(1)+atmp(3)+.015,atmp(2),.02,atmp(4)])
%[~,idl] = sort(idx);
itmp = groupColors(cIdx,:);
itmp = permute(repmat(itmp,1,1,4),[1,3,2]);
colorIDs = itmp;%(idl,:,:);
whiteStrips = zeros(size(colorIDs));
whiteStrips(1,3:4,:)=1;
stripEnd=true;
for j=2:size(whiteStrips,1) % this is a dumb way to do this
    if (colorIDs(j-1,1,1)==colorIDs(j,1,1))&&(colorIDs(j-1,1,2)==colorIDs(j,1,2))&&(colorIDs(j-1,1,3)==colorIDs(j,1,3))
        whiteStrips(j,:,:)=whiteStrips(j-1,:,:);
    else
        stripEnd=~stripEnd;
        if stripEnd; whiteStrips(j,3:4,:)=1;
        else; whiteStrips(j,1:2,:)=1;
        end
    end
end
colorIDs(find(whiteStrips))=0;
imagesc(colorIDs)
setFigColors;
axis off

f2.InvertHardcopy = 'off';
if ~isfolder([savePath,'_plots/']); mkdir([savePath,'_plots/',fnm]); end
expNameHandle=strrep(expID,'/','_');
saveas(f2, [savePath,'_plots/',expNameHandle,'_',fnm,'_from',ccFlag(end-2:end),'.png'])



function setFigColors
set(gca,'Fontsize',14)
set(gca,'color','none');%'k')
set(gcf,'color','none');%'k');
set(gca,'xcolor',[.9 .9 .9]);
set(gca,'ycolor',[.9 .9 .9]);