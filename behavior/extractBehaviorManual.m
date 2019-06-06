%baseDir = '/Volumes/RASPIUSB/092818_looming/';%'/Volumes/dataFast/habaRegistered/'; %'/Volumes/SCAPEdata1/scapeBehavior/';
%date = 'fly4_CS'; %'2018_08_24_odor'; %'blue_mp4s';

function extractBehaviorManual(date,behaviorFolder)
% loops through videos, checks ROI location, and extracts light/shock timeseries
global lsf editROI waitForPress vidNum videoFolder files tStep
global shapeInserterBubble shapeInserterLegs viewROI legs RGB J isavi
global playFobj playRobj imagObj stimObj drinkObj editObj ffRobj ffFobj
global imagIsOn stimIsOn drinkIsOn editIsOn

if ~exist('date','var')
    date = '0828_f2r3/';%'2019_02_14_Nsyb_NLS6s_Su/'; %'2018_10_24_MCFO_IRtest/';%'2018_11_1_looming/';%'2018_08_24_odorAndRun/';
    behaviorFolder = '/Users/evan/Dropbox/_sandbox/sourceExtraction/good/'; %'/Users/evan/Dropbox/_sandbox/behavTest/';%'/Volumes/dataFast2/scapeBehavior/';
end
%behaviorFolder = '/Users/evan/Documents/_Lab/Data/_Behavior/'; %'/Users/evan/Dropbox/_temp/'; %
parentFolder = [behaviorFolder,date]; %behaviorFolder; %
videoFolder = parentFolder; %[parentFolder,'UncompressedAVI/'];

isavi = false; %true;
if isavi
    files = dir([videoFolder,'f*.avi']);
else
    files = dir([videoFolder,'f*.mp4']);
end

roiTemplate = load([behaviorFolder,'_ROItemplate.mat']);


% choose ROI's and save matfile

waitForPress = true;
%for n = 1:length(avifiles)
vidNum = 1;
movfile = [videoFolder,files(vidNum).name];

imagIsOn = 0;
stimIsOn = 0;
drinkIsOn = 0;
editIsOn = 0;

%if waitForPress
tStep=101;
aviobj = VideoReader(movfile);
for t=1:100; readFrame(aviobj); end % dump first 100 frames
objectFrame = readFrame(aviobj);
RGB = objectFrame; %repmat(objectFrame,[1,1,3]); % convert I to an RGB image
lsf = figure;
set(gcf,'color','k')
pos=get(gcf,'Position');
heightRescale = 2;
widthRescale = 2;
pos(3)=pos(3)*widthRescale;
pos(4)=pos(4)*heightRescale;
set(gcf,'Position',pos)
axes('Position',[.05 .12 .8 .8])

uicontrol('Style', 'pushbutton', 'String', 'Extract',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.875 .59 .1 .1],...
    'Callback', @extractButtonCallback);        % Pushbutton string callback
uicontrol('Style', 'pushbutton', 'String', 'Fix',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.875 .85 .1 .1],...
    'Callback', @fixButtonCallback);        % Pushbutton string callback
uicontrol('Style', 'pushbutton', 'String', 'Accept',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.875 .72 .1 .1],...
    'Callback', @acceptButtonCallback);        % Pushbutton string callback
uicontrol('Style', 'pushbutton', 'String', 'SAVE',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.875 .46 .1 .1],...
    'Callback', @saveButtonCallback);        % Pushbutton string callback


imagObj = uicontrol('Style', 'togglebutton', 'String', 'imaging',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.875 .3 .1 .08],...
    'Callback', @imagingButtonCallback);        
stimObj = uicontrol('Style', 'togglebutton', 'String', 'stimulus',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.875 .2 .1 .08],...
    'Callback', @stimButtonCallback);        
drinkObj = uicontrol('Style', 'togglebutton', 'String', 'drinking',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.875 .1 .1 .08],...
    'Callback', @drinkButtonCallback);        

uicontrol('Style', 'pushbutton', 'String', '<',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.1 .11 .1 .04],...
    'Callback', @stepBackButtonCallback);        % Pushbutton string callback
uicontrol('Style', 'pushbutton', 'String', '>',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.7 .11 .1 .04],...
    'Callback', @stepForwardButtonCallback);        % Pushbutton string callback
playRobj = uicontrol('Style', 'togglebutton', 'String', '<<',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.1 .065 .1 .04],...
    'Callback', @playBackButtonCallback);        
playFobj = uicontrol('Style', 'togglebutton', 'String', '>>',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.7 .065 .1 .04],...
    'Callback', @playForwardButtonCallback);        
ffRobj = uicontrol('Style', 'togglebutton', 'String', '<<<',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.1 .02 .1 .04],...
    'Callback', @ffBackButtonCallback);        
ffFobj = uicontrol('Style', 'togglebutton', 'String', '>>>',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.7 .02 .1 .04],...
    'Callback', @ffForwardButtonCallback);        

editObj = uicontrol('Style', 'togglebutton', 'String', 'edit',...
    'Units','normalized',...
    'Fontsize',14,...
    'ForegroundColor',[0 0 0],...
    'Position', [.4 .02 .1 .1],...
    'Callback', @editButtonCallback);        


yellow = uint8([255 255 0]); % [R G B]; class of yellow must match class of I
cyan = uint8([0 255 255]); 
shapeInserterBubble = vision.ShapeInserter('BorderColor','Custom','CustomBorderColor',yellow);
shapeInserterLegs = vision.ShapeInserter('BorderColor','Custom','CustomBorderColor',cyan);
viewROI = roiTemplate.viewROI;
legs = roiTemplate.legs; %zeros(1,4);
try
    J = step(shapeInserterBubble,RGB,int32(viewROI));
catch
    J = step(shapeInserterBubble,cat(3,RGB,RGB,RGB),int32(viewROI));
end
J = step(shapeInserterLegs,J,int32(legs));
imshow(J);
title(files(vidNum).name,'fontsize',20,'color','w','interpreter','none');


editROI = true;

while waitForPress
    pause(.1)
end

while editROI
    viewROI = zeros(size(viewROI));
    legs = zeros(size(viewROI));
    try
        J = step(shapeInserterBubble,RGB,int32(viewROI));
    catch
        J = step(shapeInserterBubble,cat(3,RGB,RGB,RGB),int32(viewROI));
    end
    J = step(shapeInserterLegs,J,int32(legs));
    imshow(J);
    
    title('Draw viewer ROI','fontsize',20,'color','w');
    waitForPress = true;
    shapeInserterBubble = vision.ShapeInserter('BorderColor','Custom','CustomBorderColor',yellow);
    viewROI=round(getPosition(imrect));
    
    title('Draw Legs ROI','fontsize',20,'color','w');
    shapeInserterLegs = vision.ShapeInserter('BorderColor','Custom','CustomBorderColor',cyan);
    legs=round(getPosition(imrect));
    
    editROI = false;
end

save( [behaviorFolder,'_ROItemplate.mat'],'viewROI','legs');




% ------- BUTTONS ------------------------------------------------------

function updateObjectsAndTraces(steps)
global tStep editIsOn traces imagObj stimObj drinkObj imagIsOn stimIsOn drinkIsOn
if ~nargin; steps=1; end
if editIsOn
    % if edit is on, use obj value to update traces
    if abs(steps)>1
        if steps<0
            trange = tStep:tStep-steps-1;
        else
            trange = tStep-steps+1:tStep;
        end
        traces.isStimOn(trange) = stimIsOn;
        traces.isImagingOn(trange) = imagIsOn;
        traces.isDrinking(trange) = drinkIsOn;
    else
        traces.isStimOn(tStep) = stimIsOn;
        traces.isImagingOn(tStep) = imagIsOn;
        traces.isDrinking(tStep) = drinkIsOn;
    end
else
    % if edit is off, use traces to update obj value
    updateObjectsFromTraces(stimObj,traces.isStimOn(tStep));
    updateObjectsFromTraces(imagObj,traces.isImagingOn(tStep));
    updateObjectsFromTraces(drinkObj,traces.isDrinking(tStep));
    stimIsOn=traces.isStimOn(tStep);
    imagIsOn=traces.isImagingOn(tStep);
    drinkIsOn=traces.isDrinking(tStep);
end

function updateObjectsFromTraces(ob,trV)
if trV
    ob.ForegroundColor = [1 0 0];
    ob.FontWeight = 'bold';
else
    ob.ForegroundColor = [0 0 0];
    ob.FontWeight = 'normal';
end

function saveButtonCallback(hObj,event)
global m traces
m.traces = traces;

function imagingButtonCallback(hObj,event)
global editIsOn tStep traces imagIsOn
if editIsOn
    if hObj.Value
        hObj.ForegroundColor = [1 0 0];
        hObj.FontWeight = 'bold';
        traces.isImagingOn(tStep) = 1;
        imagIsOn = 1;
    else
        hObj.ForegroundColor = [0 0 0];
        hObj.FontWeight = 'normal';
        traces.isImagingOn(tStep) = 0;
        imagIsOn = 0;
    end
end

function stimButtonCallback(hObj,event)
global editIsOn tStep traces stimIsOn
if editIsOn
    if hObj.Value
        hObj.ForegroundColor = [1 0 0];
        hObj.FontWeight = 'bold';
        traces.isStimOn(tStep) = 1;
        stimIsOn = 1;
    else
        hObj.ForegroundColor = [0 0 0];
        hObj.FontWeight = 'normal';
        traces.isStimOn(tStep) = 0;
        stimIsOn = 0;
    end
end

function drinkButtonCallback(hObj,event)
global editIsOn tStep traces drinkIsOn
if editIsOn
    if hObj.Value
        hObj.ForegroundColor = [1 0 0];
        hObj.FontWeight = 'bold';
        traces.isDrinking(tStep) = 1;
        drinkIsOn = 1;
    else
        hObj.ForegroundColor = [0 0 0];
        hObj.FontWeight = 'normal';
        traces.isDrinking(tStep) = 0;
        drinkIsOn = 0;
    end
end

function editButtonCallback(hObj,event)
global editIsOn
if hObj.Value
    hObj.ForegroundColor = [0 .8 0];
    hObj.FontWeight = 'bold';
    editIsOn = 1;
else
    hObj.ForegroundColor = [0 0 0];
    hObj.FontWeight = 'normal';
    editIsOn = 0;
end

function stepBackButtonCallback(hObj,event)
global im tStep
[w,~,nframe]=size(im);
tStep = tStep-1;
imshow(im(:,:,tStep));
text(20,w-20,['Frame ',num2str(tStep),' of ',num2str(nframe)],'color',[1 1 1]);
updateObjectsAndTraces;

function stepForwardButtonCallback(hObj,event)
global im tStep
[w,~,nframe]=size(im);
tStep = tStep+1;
imshow(im(:,:,tStep));
text(20,w-20,['Frame ',num2str(tStep),' of ',num2str(nframe)],'color',[1 1 1]);
updateObjectsAndTraces;

function playBackButtonCallback(hObj,event)
global isActive
isActive=1;
if hObj.Value
    hObj.ForegroundColor = [0 .8 0];
    hObj.FontWeight = 'bold';
    playR;
else
    isActive=0;
    hObj.ForegroundColor = [0 0 0];
    hObj.FontWeight = 'normal';
end

function playForwardButtonCallback(hObj,event)
global isActive
isActive=1;
if hObj.Value
    hObj.ForegroundColor = [0 .8 0];
    hObj.FontWeight = 'bold';
    playF;
else
    hObj.ForegroundColor = [0 0 0];
    hObj.FontWeight = 'normal';
    isActive=0;
end

function playF
global im tStep isActive playFobj
[w,~,n]=size(im);
while isActive && tStep<n
    tStep = tStep+1;
    imshow(im(:,:,tStep));
    text(20,w-20,['Frame ',num2str(tStep),' of ',num2str(n)],'color',[1 1 1],'fontweight','bold');
    updateObjectsAndTraces;
    drawnow
end
if tStep==n
    isActive=false;
    playFobj.ForegroundColor = [0 0 0];
    playFobj.FontWeight = 'normal';
end

function playR
global im tStep isActive playRobj
[w,~,n]=size(im);
while isActive && tStep>1
    tStep = tStep-1;
    imshow(im(:,:,tStep));
    text(20,w-20,['Frame ',num2str(tStep),' of ',num2str(n)],'color',[1 1 1],'fontweight','bold');
    updateObjectsAndTraces;
    drawnow
end
if tStep==1
    isActive=false;
    playRobj.ForegroundColor = [0 0 0];
    playRobj.FontWeight = 'normal';
end

function ffBackButtonCallback(hObj,event)
global isActive
isActive=1;
if hObj.Value
    hObj.ForegroundColor = [0 .8 0];
    hObj.FontWeight = 'bold';
    ffR;
else
    isActive=0;
    hObj.ForegroundColor = [0 0 0];
    hObj.FontWeight = 'normal';
end

function ffForwardButtonCallback(hObj,event)
global isActive
isActive=1;
if hObj.Value
    hObj.ForegroundColor = [0 .8 0];
    hObj.FontWeight = 'bold';
    ffF;
else
    hObj.ForegroundColor = [0 0 0];
    hObj.FontWeight = 'normal';
    isActive=0;
end

function ffF
global im tStep isActive ffFobj
[w,~,n]=size(im);
while isActive && tStep<n
    st = min(10,n-tStep);
    tStep = tStep+st;
    imshow(im(:,:,tStep));
    text(20,w-20,['Frame ',num2str(tStep),' of ',num2str(n)],'color',[1 1 1],'fontweight','bold');
    updateObjectsAndTraces(st);
    drawnow
end
if tStep==n
    isActive=false;
    ffFobj.ForegroundColor = [0 0 0];
    ffFobj.FontWeight = 'normal';
end

function ffR
global im tStep isActive ffRobj
[w,~,n]=size(im);
while isActive && tStep>1
    st = min(10,tStep-1);
    tStep = tStep-st;
    imshow(im(:,:,tStep));
    text(20,w-20,['Frame ',num2str(tStep),' of ',num2str(n)],'color',[1 1 1],'fontweight','bold');
    updateObjectsAndTraces(-st);
    drawnow
end
if tStep==1
    isActive=false;
    ffRobj.ForegroundColor = [0 0 0];
    ffRobj.FontWeight = 'normal';
end














function extractButtonCallback(hObj,event)
global editROI waitForPress lsf tStep
global videoFolder files isavi im m traces


editROI = false;
waitForPress = false;
%close(lsf)
title('EXTRACTING','fontsize',20,'color','r');

for n = 1:length(files)
    traces = struct;
    display(['Tracking ',num2str(n),' of ', num2str(length(files))]);
    movfile = [videoFolder,files(n).name];
    m = matfile( [movfile(1:end-3),'mat'],'Writable',true);
    aviobj = VideoReader(movfile);
    nframe = floor(aviobj.Duration*aviobj.FrameRate); %aviobj.NumberOfFrames;
    rb = m.viewROI;
    rl = m.legs;

    legMean = zeros(nframe,1);
    legVar = zeros(nframe,1);
    
    legFrame = zeros(rl(4)+1,rl(3)+1,2);
    im = zeros(rb(4)+1,rb(3)+1,nframe,'uint16');
    
    if isavi
        frame = im2uint16(readFrame(aviobj));
    else
        frame = rgb2gray(im2uint16(readFrame(aviobj)));
    end
    legFrame(:,:,1) = frame(rl(2):rl(2)+rl(4),rl(1):rl(1)+rl(3));
    im(:,:,1) = frame(rb(2):rb(2)+rb(4),rb(1):rb(1)+rb(3));
    
    %nframe=200;
    %disp('only loaded 200 frames')
    for t=2:nframe
        if ~mod(t,50); disp(t/nframe); end
        if isavi
            frame = im2uint16(readFrame(aviobj));
        else
            frame = rgb2gray(im2uint16(readFrame(aviobj)));
        end
        im(:,:,t) = frame(rb(2):rb(2)+rb(4),rb(1):rb(1)+rb(3));

        legFrame(:,:,2) = frame(rl(2):rl(2)+rl(4),rl(1):rl(1)+rl(3));
        legMean(t) = mean(mean(legFrame(:,:,2))); %mean(mean(frame));
        legVar(t) = mean(mean( (legFrame(:,:,2)-legFrame(:,:,1)).^2 ))/mean(mean(legFrame(:,:,2).^2));
        legFrame(:,:,1) = legFrame(:,:,2);
    end
    traces.legMean = legMean;
    traces.legVar = legVar;
    traces.isImagingOn = zeros(size(legMean));
    traces.isStimOn = zeros(size(legMean));
    traces.isDrinking = zeros(size(legMean));
    m.traces = traces;
    
    %figure; plot(bubMean); drawnow;
    %time.imOn = input('imaging onset:');       % fix this (not usable anyway for looming dataset)
    %time.imOff = input('imaging offset:');
    %m.time = time;
    
    %m.time.vps = input('volume rate:');
    %m.time.behavTime = (1:(m.time.imOff-m.time.imOn))/m.time.vps;
end
imshow(im(:,:,tStep));
title('annotate manually','fontsize',20,'color','w');
w=size(im,1);
text(20,w-20,['Frame ',num2str(tStep),' of ',num2str(nframe)],'color',[1 1 1],'fontweight','bold');




function fixButtonCallback(hObj,event)
global editROI waitForPress

editROI = true;
waitForPress = false;



function acceptButtonCallback(hObj,event)
global vidNum videoFolder files shapeInserterBubble viewROI legs RGB J

save( [videoFolder,files(vidNum).name(1:end-3),'mat'],'viewROI','legs');
% vidNum = vidNum+1;
% if vidNum<=length(files)
%     movfile = [videoFolder,files(vidNum).name];
%     
%     aviobj = VideoReader(movfile);
%     for t=1:100; readFrame(aviobj); end  % dump first 100 frames
%     objectFrame = readFrame(aviobj);
%     RGB = objectFrame; %repmat(objectFrame,[1,1,3]); % convert I to an RGB image
%     
%     yellow = uint8([255 255 0]); % [R G B]; class of yellow must match class of I
%     shapeInserterBubble = vision.ShapeInserter('BorderColor','Custom','CustomBorderColor',yellow);
%     %viewROI = roiTemplate.viewROI;
%     try
%         J = step(shapeInserterBubble,RGB,int32(viewROI));
%     catch
%         J = step(shapeInserterBubble,cat(3,RGB,RGB,RGB),int32(viewROI));
%     end
%     J = step(shapeInserterTrigger,J,int32(trig));
%     J = step(shapeInserterLegs,J,int32(legs));
%     imshow(J);
%     title(files(vidNum).name,'fontsize',20,'color','w','interpreter','none');
% else
%     title('FILES COMPLETE','fontsize',20,'color','r');
% end

