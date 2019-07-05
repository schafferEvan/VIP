
function extractBehaviorAuto(videoFolder)
% loops through videos, checks ROI location, and extracts light/shock timeseries

addpath(genpath('..'))

% if ~exist('date','var')
%     date = '0828_f2r3/';%'2019_02_14_Nsyb_NLS6s_Su/'; %'2018_10_24_MCFO_IRtest/';%'2018_11_1_looming/';%'2018_08_24_odorAndRun/';
%     behaviorFolder = '/Users/evan/Dropbox/_sandbox/sourceExtraction/good/'; %'/Users/evan/Dropbox/_sandbox/behavTest/';%'/Volumes/dataFast2/scapeBehavior/';
% end
% parentFolder = [behaviorFolder,date]; %behaviorFolder; %
% videoFolder = parentFolder; %[parentFolder,'UncompressedAVI/'];

% parameters ------------------------------------------------------------
ballthresh = 0.6;       % pixel threshold for ball vs not ball (quantile of blurred image)
nframes = 1000;         % num frames from which to estimate ball roi
indicatorHeight = 10;   % number of rows at top of image from which to measure indicator signal
% ------------------------------------------------------------------------

avifiles = dir([videoFolder,'f*.avi']);
mp4files = dir([videoFolder,'f*.mp4']);

if ~isempty(avifiles) 
    if ~isempty(mp4files)
        error('directory contains both avi and mp4 files')
    else
        files = avifiles;
    end
else
    files = mp4files;
end

for j=1:length(files)
    movfile = [videoFolder,files(j).name];
    ballROI = getBallROI(movfile, ballthresh, nframes); % find ROI of ball
    extract(movfile, ballROI, indicatorHeight)          % extract ball motion energy and indicator signal (1st PC of top of image)
end


function ballROI = getBallROI(movfile, ballthresh, nframes)
aviobj = VideoReader(movfile);

ballAvg = zeros(aviobj.Height,aviobj.Width);
for t=1:nframes
    frame = rgb2gray(im2uint16(readFrame(aviobj))); 
    ballAvg = ballAvg+double(frame)/nframes; 
end


% find connected regions
imTmp = conv2(ballAvg,1/50^2*ones(50),'same');
imthresh = quantile(imTmp(:),ballthresh);

ImAvgbw = ballAvg>imthresh;
cc=bwconncomp(ImAvgbw,6);
blobStats=regionprops(cc,'Area','Centroid','ConvexImage','BoundingBox');

% find largest connected region (ball)
M=0;id=0; 
for j=1:length(blobStats)
    if blobStats(j).Area>M
        id=j; 
        M=blobStats(j).Area; 
    end
end

% define ballROI as the convex hull of the largest connected region
ballbw = zeros(size(ImAvgbw));
vExtent = floor(blobStats(id).BoundingBox(2))+(0:blobStats(id).BoundingBox(4)-1);
hExtent = floor(blobStats(id).BoundingBox(1))+(0:blobStats(id).BoundingBox(3)-1);
ballbw(vExtent, hExtent) = blobStats(id).ConvexImage;
ballROI = find(ballbw);




function extract(movfile, ballROI, indicatorHeight)
global isavi m traces


m = matfile( [movfile(1:end-3),'mat'],'Writable',true);
aviobj = VideoReader(movfile);
nframe = floor(aviobj.Duration*aviobj.FrameRate); %aviobj.NumberOfFrames;
m.ball = ballROI;
m.indicatorHeight = indicatorHeight;

legMean = zeros(nframe,1);
legVar = zeros(nframe,1);
indicatorMat = zeros(aviobj.Width*indicatorHeight, nframe);

legFrame = zeros(length(ballROI),2);

if isavi
    frame = im2uint16(readFrame(aviobj));
else
    frame = rgb2gray(im2uint16(readFrame(aviobj)));
end
legFrame(:,1) = frame(ballROI);
indicatorMat(:,1) = reshape( frame(1:indicatorHeight,:,1), aviobj.Width*indicatorHeight, 1 );

% to visualize ball ROI:
% f=frame;
% f(ballROI)=8e4
% figure; imagesc(f)

for t=2:nframe
    if ~mod(t,50); disp(t/nframe); end
    if isavi
        frame = im2uint16(readFrame(aviobj));
    else
        frame = rgb2gray(im2uint16(readFrame(aviobj)));
    end
    indicatorMat(:,t) = reshape( frame(1:indicatorHeight,:,1), aviobj.Width*indicatorHeight, 1 );
    
    legFrame(:,2) = frame(ballROI);
    legMean(t) = mean(mean(legFrame(:,2))); %mean(mean(frame));
    legVar(t) = mean(mean( (legFrame(:,2)-legFrame(:,1)).^2 ))/mean(mean(legFrame(:,2).^2));
    legFrame(:,1) = legFrame(:,2);
end
traces.legMean = legMean;
traces.legVar = legVar;
traces.isImagingOn = zeros(size(legMean));
traces.isStimOn = zeros(size(legMean));
traces.isDrinking = zeros(size(legMean));

%c = pca(indicatorMat);
indicatorMat=indicatorMat-mean(indicatorMat(:));
traces.isImagingOn = mean(indicatorMat,1);
m.traces = traces;

