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



% generate xticks and labels
tks = floor(time(1)/60):floor(time(end)/60);
tkLoc = zeros(size(tks));
tkLab = cell(size(tks));
labelFlagVal = round(length(tks)/4); % only label 4 ticks total
for j=1:length(tks)
    [~,tkLoc(j)] = min(abs(time-60*tks(j)));
    tkLab{j}=num2str(tks(j));
end
for j=1:length(tks)
    if mod(j-1,labelFlagVal)
        tkLab{length(tks)-j+1}='';
    end
end
%tkLab{length(tks)}=num2str(length(tks));
% tkLab{1}='';

f2 = figure;
set(gcf,'color','w')
pos=get(gcf,'Position');
pos(3)=500;%pos(3)*1.5;
pos(4)=800;%pos(4)*1.5;
set(gcf,'Position',pos)
f2.InvertHardcopy = 'off';
f2.PaperUnits = 'points';
f2.PaperSize = 1.1*[pos(3) pos(4)];

a1=subplot(5,1,2:5); imagesc(pData(pOrdFull,:)); xlim([1 size(pData,2)]); caxis([0.1 0.9])
ylabel('Cell Number','fontsize',12); xlabel('Time (min)','fontsize',12);
set(gca,'XTick',tkLoc); set(gca,'XTickLabel',tkLab);
box off
colormap hot
setFigColors;
tmpA=get(gca,'Position');

totNheaderPlots = fd.showBallVar+fd.showDrink+fd.showDLC*size(fd.dlcData,2)/3;
headerBottom = 1.05*(a1.Position(2)+a1.Position(4));
headerTotHeight = 0.92 - headerBottom;
headerFigHeight = headerTotHeight/totNheaderPlots;
nhPlotted = 0;
expTitle=strrep(expID,'_','\_');
headerLabelXoffset = -0.09*size(pData,2);

disp([num2str(totNheaderPlots),' header plots'])

if fd.showDrink
    axes('position',[a1.Position(1), headerBottom+headerFigHeight*(totNheaderPlots-1), a1.Position(3), headerFigHeight])
    plot(alignedBehavior.stim,'linewidth',1); hold all; plot(alignedBehavior.drink,'linewidth',1); 
    xlim([1 size(pData,2)]); ylim([0 1]); 
    axis off
    %tmpC=get(gca,'Position');
    %set(gca,'Position',[tmpC(1),tmpC(2)-.4*(tmpC(2)-tmpB(2)-tmpB(4)),tmpC(3),tmpC(4)+.4*(tmpC(2)-tmpB(2)-tmpB(4))]);
    nhPlotted = nhPlotted + 1;
    title([expTitle,dFlag,ccFlag],'color',[.9,.9,.9]);%,'interpreter','none')
    text(headerLabelXoffset, 0.5, 'drink', 'Fontsize',8, 'color',[.9,.9,.9])
end
if fd.showBallVar
    axes('position',[a1.Position(1), headerBottom+headerFigHeight*(totNheaderPlots-1-nhPlotted), a1.Position(3), headerFigHeight]); 
    plot(behNorm,'color',[.9,.9,.9],'linewidth',1); 
    xlim([1 size(pData,2)]); ylim([min(behNorm) max(behNorm)]); 
    axis off
    %tmpB=get(gca,'Position');
    %set(gca,'Position',[tmpB(1),tmpB(2)-.4*(tmpB(2)-tmpA(2)-tmpA(4)),tmpB(3),tmpB(4)+.4*(tmpB(2)-tmpA(2)-tmpA(4))]);
    if ~nhPlotted
        title([expTitle,dFlag,ccFlag],'color',[.9,.9,.9]);%,'interpreter','none')
    end
    nhPlotted = nhPlotted + 1;
    text(headerLabelXoffset, .5*(max(behNorm)-min(behNorm)), 'ball', 'Fontsize',8, 'color',[.9,.9,.9])
end

dlcColor = jet(size(fd.dlcData,2)/3);
for j=1:totNheaderPlots-nhPlotted
    axes('position',[a1.Position(1), headerBottom+headerFigHeight*(totNheaderPlots-1-nhPlotted), a1.Position(3), headerFigHeight]); 
    xdataChunk = diff(fd.dlcData(:,1+(j-1)*3)); %.(['x',num2str(j)]);
    ydataChunk = diff(fd.dlcData(:,2+(j-1)*3)); %.(['y',num2str(j)]);
    legEnergy = xdataChunk.^2 + ydataChunk.^2;
    legEnergy = TVL1denoise(legEnergy, 0.1, 100); 
    hold on; plot(legEnergy,'color',dlcColor(j,:),'linewidth',1);
    xlim([1 size(pData,2)]); %ylim([min(behNorm) max(behNorm)]); 
    axis off
    if ~nhPlotted
        title([expTitle,dFlag,ccFlag],'color',[.9,.9,.9]);%,'interpreter','none')
    end
    nhPlotted = nhPlotted + 1;
    text(headerLabelXoffset, nanmean(legEnergy), ['point ',num2str(j)], 'Fontsize',8, 'color',dlcColor(j,:))
end




% make colorbar to relate footprints to traces ----------------------------
cIdx = zeros(size(idx));
tmpIdx = 0;
for k=1:K
    kIds = idx==pOrd(k); %goodIds(idx==pOrd(k));
    cIdx(tmpIdx+(1:sum(kIds)))=k; tmpIdx=tmpIdx+sum(kIds);
end
groupColors = flipud(jet(K));

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
if ~isfolder([savePath,'_plots/']); mkdir([savePath,'_plots/']); end
expNameHandle=strrep(expID,'/','_');
saveas(f2, [savePath,'_plots/',expNameHandle,'_',fnm,'_from',ccFlag(end-2:end),'.png'])



function setFigColors
set(gca,'Fontsize',14)
set(gca,'color','none');%'k')
set(gcf,'color','none');%'k');
set(gca,'xcolor',[.9 .9 .9]);
set(gca,'ycolor',[.9 .9 .9]);