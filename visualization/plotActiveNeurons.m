
function plotActiveNeurons(codePath, expDir, expID, savePath, fromGreenCC, showBallVar, showDrink, showDLC)

if ~nargin
    codePath = '..';
    expDir = '/Volumes/SCAPEdata1/finalData/2019_07_01_Nsyb_NLS6s_walk/fly2/Yproj/';
    expID = '2019_07_01_Nsyb_NLS6s_walk/fly2';
    savePath = '/Volumes/SCAPEdata1/figsAndMovies/2019_07_01_Nsyb_NLS6s_walk/fly2/';
    fromGreenCC = false;    % use ROIs parsed from green image (false -> use red)
    showBallVar = true;     % show motion energy of ball extracted as pix var
    showDrink = false;      % show trace of bubble and drinking
    showDLC = true;         % show output of deeplabcut
end

if ischar(fromGreenCC); fromGreenCC = str2double(fromGreenCC); end % this can happen if called from bash script
if ischar(showBallVar); showBallVar = str2double(showBallVar); end % this can happen if called from bash script
if ischar(showDrink); showDrink = str2double(showDrink); end % this can happen if called from bash script
if ischar(showDLC); showDLC = str2double(showDLC); end % this can happen if called from bash script

fd.showBallVar = showBallVar;
fd.showDrink = showDrink;
fd.showDLC = showDLC;

addpath(genpath(codePath))

if fromGreenCC
    load([expDir,'/post_fromYcc.mat'])
    load([expDir,'/F.mat'])
    matfileRaw = matfile([expDir,'/F.mat']);
    matfilePost = matfile([expDir,'/post_fromYcc.mat']);
    ccFlag = ' from Ycc';
else
    load([expDir,'/post_fromRcc.mat'])
    load([expDir,'/F_fromRed.mat'])
    matfileRaw = matfile([expDir,'/F_fromRed.mat']);
    matfilePost = matfile([expDir,'/post_fromRcc.mat']);
    ccFlag = ' from Rcc';
end
load([expDir,'/Ysum.mat'])
load([expDir,'/alignedBehavAndStim.mat'])
load([expDir,'/alignedBehavSmooth.mat'])

fd = struct('expID',expID,'savePath',savePath,'time',time,'alignedBehavior',alignedBehavior,...
    'ccFlag',ccFlag,'A',A,'Ysum',Ysum,'showBallVar',showBallVar,'showDrink',showDrink,'showDLC',showDLC);
[fd.d1,fd.d2,fd.d3] = size(Ysum);

beh = behSmooth; %alignedBehavior.legVar; %

normalizeInPieces = false; % set to true for run/flail/run experiments
if normalizeInPieces
    behNorm = zeros(size(beh));
    trialTmp = unique(trialFlag);
    for j=1:length(trialTmp)
        behPiece = beh(trialFlag==trialTmp(j));
        behNormPiece = normalizeBehavior(behPiece);
        behNorm(trialFlag==trialTmp(j)) = behNormPiece;
    end
else
    behNorm = normalizeBehavior(beh);
end

fd.behNorm = behNorm;
fd.dlcData = alignedBehavior.dlcData;

dYYfull = dYY; %zeros(size(F));
dRRfull = dRR; %zeros(size(F));
dOOfull = dOO; %zeros(size(F));
fd.goodIds = find(matfilePost.goodIds); %ones(size(matfilePost.goodIds)); %

oData = dOOfull(fd.goodIds,:); %Fsc(goodIds,:) - F0full(goodIds,:);
yData = dYYfull(fd.goodIds,:);
rData = dRRfull(fd.goodIds,:);

fd.idx = clustInd; %(find(isGreenDominant)); %clustInd
fd.K=length(unique(fd.idx)); %size(clustData,1);
m=clustData;

m(isnan(m))=0;
maxdevs = movmad( smoothdata(m,2,'movmean',20), 100,2);

[~,seed] = max( max(maxdevs,[],2));
mc = corrcoef(m');

fd.pOrd = zeros(1,fd.K);
fd.pOrd(1)=seed;
for k=1:fd.K-1
    tmp = nanmean(mc(fd.pOrd(1:k),:),1);
    tmp(fd.pOrd(1:k))=nan;
    [~,loc] = nanmax(tmp);
    fd.pOrd(k+1)=loc;
end
pOrdFull = zeros(size(fd.idx));
ctr=0;
for k=1:fd.K
    i = find(fd.idx==fd.pOrd(k));
    pOrdFull(ctr+(1:length(i))) = i;
    ctr = ctr+length(i);
end
fd.pOrdFull=pOrdFull;
% fd.pOrd=1;
% fd.pOrdFull = ones(1,size(rData,1));

for flist = 1:3
    if flist==1
        fd.pData = oData;
        fd.dFlag = ' \Delta O/O';
        fd.fnm = 'dOO';
    elseif flist==2
        fd.pData = yData;
        fd.dFlag = ' \Delta Y/Y';
        fd.fnm = 'dYY';
    else
        fd.pData = rData;
        fd.dFlag = ' \Delta R/R';
        fd.fnm = 'dRR';
    end
    makeRasterFig(fd);
end
makeFootprintFig(fd);



function behNorm = normalizeBehavior(beh)
cmax = quantile(beh,.9); %quantile(beh,.999); %max(beh);
cmin = quantile(beh,.1); %quantile(beh,.1); %min(beh);
behNorm = (beh-cmin)/(cmax-cmin);
behNorm(behNorm>1)=1;
behNorm(behNorm<0)=0;

