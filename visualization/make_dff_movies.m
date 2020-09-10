
addpath(genpath('~/Dropbox/_code/'))
% tracesFolder = '/Volumes/dataFast2/rustyOut/2019_02_07_Nsyb_NLS6s_Su/fly2/'; %'/Volumes/dataFast2/rustyOut/2019_01_10_hunger/';%'/Volumes/dataFast/rustyOut/2018_08_24_odor_fly3/'; %'~/Dropbox/_code/_dataTemp/pythonCaiman/run1_asis/';%'~/Dropbox/_AxelLab/_data/_scape/data/fly1nucTraces/';%'~/Dropbox/_AxelLab/_data/_scape/data/test/'; %
params.trialName = '2018_08_24_NsybNLS_odors_fly3_run1';%'reg20180824NsybNLSodorsGfly3';
% backgroundFile = [tracesFolder,params.trialName,'fr20241.tiff'];%[tracesFolder,params.trialName,'fr1181.tiff'];%'/Users/evan/caiman_data/example_movies/tiffs/reg20180801IRtestfly1fr1.tiff';%'~/Desktop/reg20180801IRtestfly1fr1.tiff'; %'~/Dropbox/_AxelLab/_data/_scape/data/fly1nucTraces/';%'~/Dropbox/_AxelLab/_data/_scape/data/test/'; %
savePath = '~/Dropbox/_AxelLab/_data/_scape/movies/';
% infoFile = dir([tracesFolder,'info/fly2_run2_info.mat']);
% load([tracesFolder,'info/',infoFile.name],'info');
%corrMat = load([tracesFolder,'imBehavCorr.mat'],'c');

load('/Users/evan/Desktop/2018_08_24_NsybNLS_odors_fly3_run1.mat')
load('/Users/evan/Desktop/2018_08_24_NsybNLS_odors_fly3_run1_Agood.mat')
load('/Users/evan/Desktop/fly3_run1_info.mat')
m = matfile('/Volumes/SCAPEdata1/finalData/2018_08_24_NsybNLS_odors/fly3/run1/2018_08_24_Nsyb_NLS6s_walk_G_fly3_run1_181reg.mat');
Y = m.R(:,:,:,1);

generateParts = true;
makeCompMovie = true;
makeScaledCompMovie = false;
makeResidMovie = false;

if generateParts
    %dOO(dOO>.3)=.3;
    dOO(dOO<.05)=.05;
    log_dOO = log(dOO)-log(.05);
    log_dOO(log_dOO<.5)=0;
    dfF = log_dOO(:,1000:2500);
    
    %Y = loadtiff(backgroundFile);
    [d1,d2,d3] = size(Y);
    
    T = size(dfF,2);
    
    %cy = var(YrA,0,2);
    %cv = var(C,0,2);
    %     sA = sum(A>0);
    %     mAC = max(A,[],1)'.*C(:,1);
    %     includeIdx = find(max(dfF,[],2)<2e4);%find( (isok).*(cv>cy).*full(sA<300)'.*full(mAC<5e2).*(max( C,[],2)<1e4) );
    %excludeList = find(mean(dfF,2)>500); 
    %includeIdx(excludeList)=0;
    %includeIdx = ones(size(dfF,1),1); includeIdx(excludeList)=0; includeIdx = find(includeIdx);
    %     dfF = dfF( includeIdx, : );
    %     A = A( :, includeIdx );
    %     dfF = smoothdata(dfF,2,'movmedian',10);
    %     m = min(dfF,[],2);
    %     dfF = dfF-m; % redefine baseline as zero
    
    % ---- STIM ----------
    %load([tracesFolder,'alignedBehavAndStim.mat'],'alignedStim');
    %stimTot = alignedStim;
    stimTot = [inf,inf];
    % --------------------
    
    
    % show movie of projections
    params.savePath = savePath;
    params.options.d1 = d1;
    params.options.d2 = d2;
    params.options.d3 = d3;
    params.plotResid = false;
    params.concatenate = false;         % use this to chain together csvs
    %A = A/diag(max(A,[],1));               % normalize spatial footprints
    params.maxArray = max( dfF,[],2); disp('max calculated')
    params.minArray = min( dfF,[],2); disp('min calculated')
    %params.maxArray = max( diag(full(sum(A,1))) * dfF,[],2); disp('max calculated')
    %params.minArray = min( diag(full(sum(A,1))) * dfF,[],2); disp('min calculated')
    %params.maxArray = max(sparse(A)*sparse(dfF),[],2); disp('max calculated')
    %params.minArray = min(sparse(A)*sparse(dfF),[],2); disp('min calculated')
    params.maxY = max(params.maxArray); %2500; disp('hardcoded max'); %max(params.maxArray);
    params.minY = min(params.minArray);
    %Y = double(Y); %/params.maxY;
end
params.Ttot = 0;
params.acqRate = 8; % volumes per second

% reshape A from python ND indexing to matlab.  loop needed to reduce RAM
% Ap = sparse( d1*d2*d3, size(A,2) );
% for j=1:size(A,2)
%     if ~mod(j,100); disp(j/size(A,2)); end
%     rois = reshape(full(A(:,j)),d3,d1,d2);
%     roisR = permute(rois,[2,3,1]);
%     Ap(:,j) = sparse(reshape(roisR,d1*d2*d3,1));
% end
% if size(b,2)>1; b = sum(b,2); end
b = zeros(size(A,1),1);
f = zeros(1,size(dfF,2));
try
    params.xum = d1*info.GUIcalFactors.x_umPerPix;
catch
    %warning('no x_um field available')
    params.xum = d1*info.GUIcalFactors.xK_umPerVolt*info.daq.scanAngle/(info.daq.pixelsPerLine-1);
end
params.yum = d2*info.GUIcalFactors.y_umPerPix;
params.zum = d3*info.GUIcalFactors.z_umPerPix;

if makeCompMovie
    params.savePath = savePath;
    scale_rois = true;
    %plot4Dproj_ES(Y, [], [d1,d2,d3],params,stimTot, A, dfF, 0*b, 0*f);
    plot4Dproj_ES(Y, [], [d1,d2,d3],params,stimTot, A, dfF, scale_rois);
    params.Ttot = params.Ttot + .1*round(10*size(dfF,2)/params.acqRate);
end

