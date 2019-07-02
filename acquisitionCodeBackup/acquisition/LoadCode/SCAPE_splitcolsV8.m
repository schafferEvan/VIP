function varargout = SCAPE_splitcolsV8(varargin)
% SCAPE_SPLITCOLSV8 MATLAB code for SCAPE_splitcolsV8.fig
%      SCAPE_SPLITCOLSV8, by itself, creates a new SCAPE_SPLITCOLSV8 or
%      raises the existingparfor
%      singleton*.
%
%      H = SCAPE_SPLITCOLSV8 returns the handle to a new SCAPE_SPLITCOLSV8 or the handle to
%      the existing singleton*.
%
%      SCAPE_SPLITCOLSV8('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCAPE_SPLITCOLSV8.M with the given input arguments.
%
%      SCAPE_SPLITCOLSV8('Property','Value',...) creates a new SCAPE_SPLITCOLSV8 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SCAPE_splitcolsV8_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SCAPE_splitcolsV8_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above scalered3 to modify the response to help SCAPE_splitcolsV8

% Last Modified by GUIDE v2.5 02-Mar-2018 17:38:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SCAPE_splitcolsV8_OpeningFcn, ...
    'gui_OutputFcn',  @SCAPE_splitcolsV8_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before SCAPE_splitcolsV8 is made visible.
function SCAPE_splitcolsV8_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SCAPE_splitcolsV8 (see VARARGIN)

% Choose default command line output for SCAPE_splitcolsV8
handles.readgui = cell2mat(varargin(1));
handles.output = hObject;


% Load data - KP edit added 7/29/2017
saw = get(handles.readgui.saw,'Value');
dsf = str2num(get(handles.readgui.dsf,'String'));
numSecondsToLoad = get(handles.readgui.secs,'String');
if strcmp(numSecondsToLoad,'all')
    numSecondsToLoad = handles.readgui.info.info.daq.scanLength;
else
    numSecondsToLoad = str2num(numSecondsToLoad);
end

% Scale data set size by spatial bin factor
scanRate = handles.readgui.info.info.daq.scanRate;               % Volumetric Scan Rate (VPS)
numFrames = handles.readgui.info.info.daq.pixelrate*handles.readgui.info.info.daq.scanLength;% length(info.daq.scanWaveform);  % Total number of frames acquired
numScans = handles.readgui.info.info.daq.numberOfScans;          % Total number of volumes acquired
framesPerScan = handles.readgui.info.info.daq.pixelsPerLine;%floor(info.camera.framerate/info.daq.scanRate);        % Number of frames per volume


%kpedit - if stage scanning instead of glavo scanning, framesPerScan=2, so
%reset number of volumes to be 1 and for all frames to be read in
%sequence
if framesPerScan == 2
    framesPerScan = round(numFrames);
    numScans = round(numFrames);
    scanRate = 1/handles.readgui.info.info.daq.scanLength;
end


numFramesPerSpool = handles.readgui.info.zyla.numFramesPerSpool;
numDepths = handles.readgui.info.zyla.numDepths;
numLatPix = handles.readgui.info.zyla.numLatPix;


% numSpoolfilesPerVolume = ceil(handles.readgui.info.info.daq.pixelsPerLine/handles.readgui.info.zyla.numFramesPerSpool);
numSpoolfilesPerVolume = ceil(framesPerScan/numFramesPerSpool);
numSpoolfilesPerVolume1 = (framesPerScan/numFramesPerSpool);

numActual_spoolfiles = length(dir([handles.readgui.directory,'/',handles.readgui.scanName,'/*.dat']));
if numSpoolfilesPerVolume1*scanRate*numSecondsToLoad>numActual_spoolfiles
    numSecondsToLoad = floor(numActual_spoolfiles/(numSpoolfilesPerVolume1*scanRate));
    disp(sprintf('%s: Not enough spool files - possible failed run, loading max available (%d secs)', handles.readgui.scanName, numSecondsToLoad));
end

if numSpoolfilesPerVolume*scanRate*numSecondsToLoad>5000
    % Creates a list of spool file names that exist for a given run this is important as the camera generates a list that appears unintelligible
    clear namesOut tempName
    spoolCounter=0;
    for i = 1:numSpoolfilesPerVolume*scanRate*numSecondsToLoad
        spoolCounter=spoolCounter+1;
        temp = i;
        for j = 1:10
            a(i, j) = mod(temp, 10^j)/(10^(j-1));
            temp = temp-mod(temp, 10^j);
        end
        tempName = mat2str(a(i, :));
        tempName = tempName(2:end-1);
        tempName = tempName(tempName ~= ' ');
        tempName = [tempName 'spool.dat'];
        namesOut{spoolCounter} = tempName;
    end
    namesOut = [{'0000000000spool.dat'} namesOut];
    handles.readgui.filesToLoad = namesOut;
    clear namesOut tempName temp a spoolCounter
end

% add extra columns for 2 buffer rows
if (2 == handles.readgui.info.info.camera.binFactor)
    numColumns = numDepths + 1;
    numRows = handles.readgui.info.zyla.ImageSize / 2 / numColumns;
    if(mod(numRows, 1) ~= 0)
        numColumns = numDepths+2;
        numRows = handles.readgui.info.zyla.ImageSize/2/numColumns;
    end
else
    numColumns = numDepths + 2;
    numRows = numLatPix;
end

Ycrop = [1 numRows-10];
xx = [1 handles.readgui.info.info.daq.pixelsPerLine];

handles.readgui.info.info.Ycrop = Ycrop;

% pre-loaded spool files (add check that there are enough)
filesToLoad = handles.readgui.filesToLoad;

% Figures out index for first frame in each volume
frame_starts = 1:framesPerScan:framesPerScan*numScans;

% Figures out which volumes to load
if (length(numSecondsToLoad) == 1)
    volumesToLoad = 1:dsf:round(numSecondsToLoad*scanRate);
else
    volumesToLoad = numSecondsToLoad(1)*scanRate:dsf:numSecondsToLoad(2)*scanRate;
end
if isempty(volumesToLoad)==0
    
    % Figures out first frame in each volume that we PLAN TO LOAD
    startFrameInVolumeToLoad = frame_starts(volumesToLoad);
    hh2 = waitbar(1,'Working');
    disp('Pre-allocating memory, please wait ...')
    
    % Creating the index of every single frames.
    numVolumeToLoad = length(startFrameInVolumeToLoad);
    endFrameInVolumeToLoad = startFrameInVolumeToLoad+framesPerScan-3;
    
    index = zeros(framesPerScan-2,numVolumeToLoad);
    
    for i = 1:numVolumeToLoad
        index(:,i) = startFrameInVolumeToLoad(i):1:endFrameInVolumeToLoad(i);
    end
    index = reshape(index,[1,((framesPerScan-2)*numVolumeToLoad)]);
    
    %         SpoolToLoad = unique([ceil(startFrameInVolumeToLoad/numFramesPerSpool),ceil(endFrameInVolumeToLoad/numFramesPerSpool)]);
    SpoolToLoad = unique(ceil(index/numFramesPerSpool));
    numFileToLoad = length(SpoolToLoad);
    numFrameToLoad = numVolumeToLoad*framesPerScan;
    
    numPixelsToReshape = numRows * numColumns * numFramesPerSpool;
    
    startSpoolInVolumeToLoad = ceil(startFrameInVolumeToLoad/numFramesPerSpool);
    startFrameInSpoolToLoad = mod(startFrameInVolumeToLoad,numFramesPerSpool);
    
    % Load all data into this variable (lateral, depth, scan, time)
    clear SCAPE_data
    
    if get(handles.readgui.makeDouble_checkbox, 'Value')
        if volumesToLoad == 1
            SCAPE_data = zeros([Ycrop(2)-Ycrop(1)]+1, numColumns-3, numFramesPerSpool,  numFileToLoad, 'double');
        else
            SCAPE_data = zeros([Ycrop(2)-Ycrop(1)]+1, numColumns-3, numFramesPerSpool,  numFileToLoad, 'double');
        end
    else
        
        if volumesToLoad == 1
            SCAPE_data = zeros([Ycrop(2)-Ycrop(1)]+1, numColumns-3, numFramesPerSpool,  numFileToLoad, 'uint16');
        else
            SCAPE_data = zeros([Ycrop(2)-Ycrop(1)]+1, numColumns-3, numFramesPerSpool,  numFileToLoad, 'uint16');
        end
    end
    %         SCAPEDATA = zeros(numRows, numColumns-2, framesPerScan-2, numVolumeToLoad);
    
    frameLoadInSpool = zeros(numFramesPerSpool,numFileToLoad);
    for i = 1:numFileToLoad
        frameLoadInSpool(:,i) = ((SpoolToLoad(i)-1)*numFramesPerSpool+1):((SpoolToLoad(i)-1)*numFramesPerSpool+numFramesPerSpool);
    end
    frameLoadInSpool = reshape(frameLoadInSpool,[1,numFramesPerSpool*numFileToLoad]);
    index = find(ismember(frameLoadInSpool,index));
    
    disp('done with memory allocation')
    
    tic
    if volumesToLoad == 1
        for spoolFileCounter = 1:numFileToLoad
            loadfile = filesToLoad{SpoolToLoad(spoolFileCounter)};
            filePath = fullfile(handles.readgui.directory, handles.readgui.scanName, loadfile);
            %The Frame index of every frames in the spoolfile
            FrameToStart = (spoolFileCounter-1)*numFramesPerSpool;
            FrameInThisSpool = FrameToStart:(FrameToStart+numFramesPerSpool-1);
            
            FID = fopen(filePath, 'r');
            % Read all data within spool file as a 1D stream of
            % integers and add it to carrydata (usually an empty
            % matrix)
            rawData = fread(FID, 'uint16=>uint16');
            
            % close file
            fclose(FID);
            rawData = reshape(rawData(1:numPixelsToReshape), ...
                [numRows,numColumns,numFramesPerSpool]);
            if get(handles.readgui.makeDouble_checkbox,'Value')
                SCAPE_data(:,:,:,spoolFileCounter) = double(rawData(Ycrop(1):Ycrop(2),1:end-3,:));
            else
                SCAPE_data(:,:,:,spoolFileCounter) = rawData(Ycrop(1):Ycrop(2),1:end-3,:);
            end
            
        end
    else
        
        parfor_progress(numFileToLoad);
        parfor (spoolFileCounter = 1:numFileToLoad, 2)
            loadfile = filesToLoad{SpoolToLoad(spoolFileCounter)};
            filePath = fullfile(handles.readgui.directory, handles.readgui.scanName, loadfile);
            %The Frame index of every frames in the spoolfile
            FrameToStart = (spoolFileCounter-1)*numFramesPerSpool;
            FrameInThisSpool = FrameToStart:(FrameToStart+numFramesPerSpool-1);
            
            FID = fopen(filePath, 'r');
            % Read all data within spool file as a 1D stream of
            % integers and add it to carrydata (usually an empty
            % matrix)
            rawData = fread(FID, 'uint16=>uint16');
            
            % close file
            fclose(FID);
            rawData = reshape(rawData(1:numPixelsToReshape), ...
                [numRows,numColumns,numFramesPerSpool]);
            if get(handles.readgui.makeDouble_checkbox,'Value')
                SCAPE_data(:,:,:,spoolFileCounter) = double(rawData(Ycrop(1):Ycrop(2),1:end-3,:));
            else
                SCAPE_data(:,:,:,spoolFileCounter) = rawData(Ycrop(1):Ycrop(2),1:end-3,:);
            end
            parfor_progress;
        end
        parfor_progress(0);
    end
    
    if volumesToLoad == 1
        SCAPE_data = reshape(SCAPE_data,[[Ycrop(2)-Ycrop(1)]+1, numColumns-3, numFramesPerSpool*numFileToLoad]);
        SCAPE_data = reshape(SCAPE_data(:,:,index),[[Ycrop(2)-Ycrop(1)]+1, numColumns-3, framesPerScan-2, numVolumeToLoad]);
    else
        SCAPE_data = reshape(SCAPE_data,[[Ycrop(2)-Ycrop(1)]+1, numColumns-3, numFramesPerSpool*numFileToLoad]);
        SCAPE_data = reshape(SCAPE_data(:,:,index),[[Ycrop(2)-Ycrop(1)]+1, numColumns-3, framesPerScan-2, numVolumeToLoad]);
    end
    
    toc
    % load waitbar
    if strcmp(handles.readgui.pf,'pf') && strcmp(get(handles.readgui.secs,'String'),'all')
        yy = [sprintf('(pf: %d secs)',numSecondsToLoad)];
    else
        yy = '';
    end
    pause(2);
    
    % Flips even frames to account for bidirectional scan. Also waits until second volume is loaded and flipped to introduce the pixel shift.
    if (saw == 0 && (round(dsf/2)== dsf/2)==0)
        SCAPE_data(:, :, :, 1:2:end) = flipdim(SCAPE_data(:, :, :, 1:2:end), 3);
    end
    handles.inputdata = SCAPE_data;
    clear SCAPE_data
    
    try
        close(hh2)
    catch
    end
    
    infob{1} = ['Directory:' handles.readgui.directory];
    infob{2} = ['Scan name:' handles.readgui.scanName];
    ss = size(handles.inputdata);
    infob{3} = sprintf('Data size: %d (y) %d (z) %d (x)',ss(1), ss(2), ss(3));
    if length(ss)==3
        ss(4) = 1;
    end
    infob{4} = sprintf('# of timepoints: %d',ss(4));
    handles.ss=ss;
    handles.writetodesktop=0;
    
    set(handles.infobox,'String',infob);
    
    handles.viewstring = {'MIP x-y','MIP y-z'};
    set(handles.MIPview,'String',handles.viewstring);
    handles.justone=1;
    
    
    % initial values for rotation and scaling
    handles.scalepastval = str2num(get(handles.ScaleG,'String'));
    handles.rotpastval = str2num(get(handles.RotateGreen,'String'));
    
    FilePath = fullfile(handles.readgui.directory, 'RGBtransform.mat');
    
    if (2 == exist(FilePath, 'file'))
        
        load(FilePath)
        set(handles.scalered,'String',num2str(transforms.scr));
        set(handles.scalegreen,'String',num2str(transforms.scg));
        set(handles.minr,'String',num2str(transforms.minr));
        set(handles.ming,'String',num2str(transforms.ming));
        handles.proc.scalered = transforms.scr;
        handles.proc.scalegreen = transforms.scg;
        handles.proc.autosc_scalered = transforms.scr;
        handles.proc.autosc_scalegreen = transforms.scg;
        handles.proc.autosc_minred = transforms.minr;
        handles.proc.autosc_mingreen = transforms.ming;
        handles.proc.autosc_minred = transforms.minr;
        handles.proc.autosc_mingreen = transforms.ming;
        
        handles.proc.rgb.xg2 = transforms.xg;
        handles.proc.rgb.xr2 = transforms.xr;
        handles.proc.rgb.yr2 = transforms.yr;
        handles.proc.rgb.yg2 = transforms.yg;
        handles.proc.rgb.zr2 = transforms.zr;
        handles.proc.rgb.zg2 = transforms.zg;
        testtime = 1;
        clear red1 green1
        set(handles.startframe,'String',num2str(testtime));
        guidata(hObject,handles)
        
        [red1, green1] = go_Callback(hObject,handles,testtime);
        
        %initial values for intensity scaling
        %Initial intensity scaling
        scr = 1/double(max(max(max(red1))));
        scg = 1/double(max(max(max(green1))));
        minr = double(min(min(min(min(red1)))));
        ming = double(min(min(min(min(green1)))));
        scr =1/(1/scr-minr);
        scg =1/(1/scg-minr);
        
        handles.proc.scalered = scr;
        handles.proc.scalegreen = scg;
        handles.proc.autosc_scalered = scr;
        handles.proc.autosc_scalegreen = scg;
        handles.proc.autosc_minred = minr;
        handles.proc.autosc_mingreen = ming;
        handles.proc.autosc_minred = minr;
        handles.proc.autosc_mingreen = ming;
        
        set(handles.scalered,'String',num2str(scr));
        set(handles.scalegreen,'String',num2str(scg));
        set(handles.minr,'String',num2str(minr));
        set(handles.ming,'String',num2str(ming));
        
        rgb(:,:,1) = squeeze(max(scr*double(red1(:,:,3:end-3,1)-minr),[],2))';
        rgb(:,:,2) = squeeze(max(scg*double(green1(:,:,3:end-3,1)-ming),[],2))';
        rgb(:,:,3) = zeros(size(rgb(:,:,1)));
        handles.RGBMIPS = rgb;
        
        set(handles.startframe,'String',num2str(testtime));
        % Update handles structure
        handles.writetodesktop=1;
        handles.justone=0;
        
    else
        figure(100)
        clf
        testtime = round(ss(4)/2);
        temp = squeeze(max(handles.inputdata(:,:,:,testtime),[],3));
        %  temp(1:floor(end/2),:) = 1*temp(1:floor(end/2),:);
        imagesc(log(double(temp))');%imagesc(squeeze(max(moviesub_sh(1:end,:,:),[],1))');
        title('select a region on each side that is the same, left to right')
        colormap gray
        [xg, zg] = ginput_col(1);
        hold on
        plot(xg,zg,'.r');
        [xr, zr] = ginput_col(1);
        hold on
        plot(xr,zr,'.g');
        
        clf
        temp= squeeze(max(handles.inputdata(:,:,:,testtime),[],2));
        %  temp(1:floor(end/2),:) = 5*temp(1:floor(end/2),:);
        % temp(floor(end/2):end,:) = 5*temp(floor(end/2):end,:);
        imagesc(log(double(temp))');%imagesc(squeeze(max(moviesub_sh(1:end,:,:),[],1))');
        title('select a region on each side that is the same, left to right')
        colormap gray
        [xg, yg] = ginput_col(1);
        hold on
        plot(xg,yg,'.r');
        [xr, yr] = ginput_col(1);
        hold on
        plot(xr,yr,'.g');
        close
        yg = yr; % this shouldn't change between color channels unless using different lasers
        handles.proc.rgb.xg2 = xg;
        handles.proc.rgb.xr2 = xr;
        handles.proc.rgb.yr2 = yr;
        handles.proc.rgb.yg2 = yg;
        handles.proc.rgb.zr2 = zr;
        handles.proc.rgb.zg2 = zg;
        clear red1 green1
        set(handles.startframe,'String',num2str(testtime));
        guidata(hObject,handles)
        
        [red1, green1] = go_Callback(hObject,handles,testtime);
        
        %Initial intensity scaling
        scr = 1/double(max(max(max(red1))));
        scg = 1/double(max(max(max(green1))));
        minr = double(min(min(min(min(red1)))));
        ming = double(min(min(min(min(green1)))));
        scr =1/(1/scr-minr);
        scg =1/(1/scg-minr);
        
        handles.proc.scalered = scr;
        handles.proc.scalegreen = scg;
        handles.proc.autosc_scalered = scr;
        handles.proc.autosc_scalegreen = scg;
        handles.proc.autosc_minred = minr;
        handles.proc.autosc_mingreen = ming;
        handles.proc.autosc_minred = minr;
        handles.proc.autosc_mingreen = ming;
        
        set(handles.scalered,'String',num2str(scr));
        set(handles.scalegreen,'String',num2str(scg));
        set(handles.minr,'String',num2str(minr));
        set(handles.ming,'String',num2str(ming));
        
        rgb(:,:,1) = squeeze(max(scr*double(red1(:,:,3:end-3,1)-minr),[],2))';
        rgb(:,:,2) = squeeze(max(scg*double(green1(:,:,3:end-3,1)-ming),[],2))';
        rgb(:,:,3) = zeros(size(rgb(:,:,1)));
        handles.RGBMIPS = rgb;
        
        set(handles.startframe,'String',num2str(testtime));
        % Update handles structure
        handles.writetodesktop=1;
        handles.justone=0;
    end
    
    axis(handles.MainFigureSC)
    imagesc(squeeze(uint8((256)*handles.RGBMIPS)));
    setappdata(handles.MainFigureSC,'fig',handles.RGBMIPS);
    axis off
else
    disp(sprintf('NOTHING TO LOAD: %s',handles.readgui.info.info.scanName));
    keyboard
end

pause(1)
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = SCAPE_splitcolsV8_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axis(handles.MainFigureSC)
setappdata(handles.MainFigureSC,'fig',handles.RGBMIPS);
imagesc(squeeze(uint8((256)*handles.RGBMIPS)));
axis off
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in go.
function [red, green] = go_Callback(hObject,handles,timepoint)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Write to tif
% profile on

xg = handles.proc.rgb.xg2;
xr = handles.proc.rgb.xr2;
yr = handles.proc.rgb.yr2;
yg = handles.proc.rgb.yg2;
zr = handles.proc.rgb.zr2;
zg = handles.proc.rgb.zg2;

clear red green
left = floor(min([xg xr]))-1;
right = floor(size(handles.inputdata,1)-xr);

top = floor(min([zr zg]))-1;
bot = floor(size(handles.inputdata,2)-max([zr zg]));
topy = floor(min([yr yg]))-1;
boty = floor(size(handles.inputdata,3)-max([yr yg]));

if length(size(handles.inputdata))>3
    green = (handles.inputdata(round(xg+[-left:right]),round(zg+[-top:bot]),round(yg+[-topy:boty]),timepoint));
    red = (handles.inputdata(round(xr+[-left:right]),round(zr+[-top:bot]),round(yr+[-topy:boty]),timepoint));
else
    green = (handles.inputdata(round(xg+[-left:right]),round(zg+[-top:bot]),round(yg+[-topy:boty])));
    red = (handles.inputdata(round(xr+[-left:right]),round(zr+[-top:bot]),round(yr+[-topy:boty])));
end

%scale gradient because imagesplitter leads to shading in opposite
%sides of red and green images
if (get(handles.gradient,'Value') == 1)
    [grad, ~, ~] = meshgrid(1:size(red,1), 1:size(red,2), 1:size(red,3));
    grad = ((grad-1)./(size(red,1)-1))*0.6+1;
    grad = flipud(permute(grad, [2 1 3]));
    clear Y Z
    green = uint16((double(green)-100).*(grad)+100);
    red = uint16((double(red)-100).*flipud(grad)+100);
end



% rotation correction for 1st gen imagesplitter
if (str2num(get(handles.RotateGreen,'String'))~=0)
    if get(handles.MIPview,'Value')== 2
        rotval = str2num(get(handles.RotateGreen,'String'));
        green = imrotate(green,rotval,'bicubic','Crop');
    end
end

% scale correction for 1st gen imagesplitter
if (str2num(get(handles.ScaleG,'String'))~=1)
    scaleval = str2num(get(handles.ScaleG,'String'));
    green = imresize(green,scaleval*[size(green,1) size(green,2)]);
    if str2num(get(handles.ScaleG,'String')) <1
        temp = zeros(size(red));
        dif1 = size(temp,1) - size(green,1);
        dif2 =  size(temp,2) - size(green,2);
        x1 = floor(dif1/2)+1 ;
        x2 = ceil(dif1/2);
        y1 = floor(dif2/2)+1;
        y2 = ceil(dif2/2);
        temp(x1:end-x2,y1:end-y2,:)=green;
        green = temp;
    else
        dif1 = size(green,1) - size(red,1);
        dif2 =  size(green,2) - size(red,2);
        x1 = floor(dif1/2)+1; x2 = ceil(dif1/2);
        y1 = floor(dif2/2)+1; y2 = ceil(dif2/2);
        green = green(x1:end-x2,y1:end-y2,:);
    end
end

guidata(hObject,handles);


% --- Executes on button press in check.
function check_Callback(hObject, handles)
% hObject    handle to check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ss = size(handles.inputdata);

if length(handles.ss) == 3; handles.ss(4) = 1; end
startframe = str2num(get(handles.startframe,'String'));
testtime = startframe;
if handles.justone==1; testtime = max(testtime); end
guidata(hObject,handles)

handles.writetodesktop=0;
clear red1 green1

[red1, green1] = go_Callback(hObject,handles,testtime);

minr = str2num(get(handles.minr,'String'));%,'%.1f'));
ming = str2num(get(handles.ming,'String'));%,'%.1f'));
scr = str2num(get(handles.scalered,'String'));%,'%.1f'));
scg = str2num(get(handles.scalegreen(1,1),'String'));%,'%.1f'));
handles.proc.scalered = scr;
handles.proc.scalegreen = scg;
handles.proc.minred = minr;
handles.proc.mingreen = ming;

if get(handles.MIPview,'Value')== 2
    clear rgb
    rgb(:,:,1) = squeeze(max(scr*(double(red1(:,:,round(end/4):round(3*end/4),1))-minr),[],3)');
    rgb(:,:,2) = squeeze(max(scg*(double(green1(:,:,round(end/4):round(3*end/4),1))-ming),[],3)');
    rgb(:,:,3) = zeros(size(rgb(:,:,1)));
    axis(handles.MainFigureSC)
    handles.RGBMIPS = rgb;
    setappdata(handles.MainFigureSC,'fig',handles.RGBMIPS);
    imagesc(squeeze(uint8((256)*handles.RGBMIPS)));
    axis off
end
if get(handles.MIPview,'Value')== 1
    clear rgb
    rgb(:,:,1) = squeeze(max(scr*(double(red1(:,round(end/4):round(3*end/4),:,1))-minr),[],2))';
    rgb(:,:,2) = squeeze(max(scg*(double(green1(:,round(end/4):round(3*end/4),:,1))-ming),[],2))';
    rgb(:,:,3) = zeros(size(rgb(:,:,1)));
    axis(handles.MainFigureSC)
    handles.RGBMIPS = rgb;
    setappdata(handles.MainFigureSC,'fig',handles.RGBMIPS);
    imagesc(squeeze(uint8((256)*handles.RGBMIPS)));
    axis off
end

title(['frame: ',num2str(testtime)])
handles.justone = 0;
handles.writetodesktop=1;

guidata(hObject, handles)


function startframe_Callback(hObject, ~, handles)
% hObject    handle to startframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startframe as scalered3
%        str2double(get(hObject,'String')) returns contents of startframe as a double
if str2num(get(handles.startframe, 'String')) > handles.ss(4)
    set(handles.startframe, 'String',num2str(handles.ss(4)))
elseif str2num(get(handles.startframe, 'String')) < 1
    set(handles.startframe, 'String','1')
end

check_Callback(hObject, handles);


% --- Executes during object creation, after setting all properties.
function startframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MIPview.
function MIPview_Callback(hObject, eventdata, handles)
% hObject    handle to MIPview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIPview contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIPview
check_Callback(hObject, handles);


% --- Executes during object creation, after setting all properties.
function MIPview_CreateFcn(hObject, ~, handles)
% hObject    handle to MIPview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function interpval_Callback(~, ~, ~)
% hObject    handle to interpval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of interpval as scalered3
%        str2double(get(hObject,'String')) returns contents of interpval as a double


% --- Executes during object creation, after setting all properties.
function interpval_CreateFcn(hObject, ~, handles)
% hObject    handle to interpval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in interpyes.
function interpyes_Callback(~, ~, ~)
% hObject    handle to interpyes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of interpyes


% --- Executes on button press in plusright.
function plusright_Callback(hObject, eventdata, handles)
% hObject    handle to plusright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if get(handles.MIPview,'Value')==1
handles.proc.rgb.xg2 = -1+ handles.proc.rgb.xg2; % = xg;
% end
% if get(handles.MIPview,'Value')==1
% handles.proc.rgb.xg2 = 1+ handles.proc.rgb.xg2;% = xg;
% end
handles.RGBMIPS = getappdata(handles.MainFigureSC,'fig');
axis(handles.MainFigureSC);
clear rgb
rgb(:,:,1) = padarray(handles.RGBMIPS(2:end-1,2:end-1,1), [1,1]);
rgb(:,:,2) = padarray(handles.RGBMIPS(2:end-1,1:end-2,2), [1,1]);
rgb(:,:,3) = zeros(size(rgb(:,:,1)));
imagesc(squeeze(uint8((256)*rgb)));
axis off
handles.RGBMIPS = rgb;
setappdata(handles.MainFigureSC,'fig',handles.RGBMIPS);
guidata(hObject,handles)
disp('image adjusted')
%check_Callback(hObject, handles);

% --- Executes on button press in minusleft.
function minusleft_Callback(hObject, eventdata, handles)
% hObject    handle to minusleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if get(handles.MIPview,'Value')==1

handles.RGBMIPS = getappdata(handles.MainFigureSC,'fig');
handles.proc.rgb.xg2 = 1+ handles.proc.rgb.xg2; % = xg;
% end
% if get(handles.MIPview,'Value')==1
% handles.proc.rgb.xg2 = 1+ handles.proc.rgb.xg2;% = xg;
% endp = str2num(get(handles.startframe,'String'));
% set(handles.startframe,'String',num2str(p-1));

%check_Callback(hObject,handles);
axis(handles.MainFigureSC);
clear rgb
rgb(:,:,1) = padarray(handles.RGBMIPS(2:end-1,2:end-1,1),[1,1]);
rgb(:,:,2) = padarray(handles.RGBMIPS(2:end-1,3:end,2),[1,1]);
rgb(:,:,3) = zeros(size(rgb(:,:,1)));
imagesc(squeeze(uint8((256)*rgb)));
axis off
handles.RGBMIPS = rgb;
setappdata(handles.MainFigureSC,'fig',handles.RGBMIPS);
guidata(hObject,handles)
display('image adjusted')


% --- Executes on button press in plustop.
function plustop_Callback(hObject, eventdata, handles)
% hObject    handle to plustop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.RGBMIPS = getappdata(handles.MainFigureSC,'fig');
if get(handles.MIPview,'Value')==1
    handles.proc.rgb.yg2 = 1+ handles.proc.rgb.yg2;% = xg;
end
if get(handles.MIPview,'Value')==2
    handles.proc.rgb.zg2 = 1+ handles.proc.rgb.zg2;% = xg;
end
axis(handles.MainFigureSC);
clear rgb
rgb(:,:,1) = padarray(handles.RGBMIPS(2:end-1,2:end-1,1),[1,1]);
rgb(:,:,2) = padarray(handles.RGBMIPS(3:end,2:end-1,2),[1,1]);
rgb(:,:,3) = zeros(size(rgb(:,:,1)));
imagesc(squeeze(uint8((256)*rgb)));
axis off
handles.RGBMIPS = rgb;

setappdata(handles.MainFigureSC,'fig',handles.RGBMIPS);
guidata(hObject,handles)
display('image adjusted')

% --- Executes on button press in minusbottom.
function minusbottom_Callback(hObject, eventdata, handles)
% hObject    handle to minusbottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.RGBMIPS = getappdata(handles.MainFigureSC,'fig');
if get(handles.MIPview,'Value')==1
    handles.proc.rgb.yg2 = -1+ handles.proc.rgb.yg2;% = xg;
end
if get(handles.MIPview,'Value')==2
    handles.proc.rgb.zg2 = -1+ handles.proc.rgb.zg2;% = xg;
end

clear rgb
rgb(:,:,1) = padarray(handles.RGBMIPS(2:end-1,2:end-1,1),[1,1]);
rgb(:,:,2) = padarray(handles.RGBMIPS(1:end-2,2:end-1,2),[1,1]);
rgb(:,:,3) = zeros(size(rgb(:,:,1)));
axis(handles.MainFigureSC);
imagesc(squeeze(uint8((256)*rgb)));
axis off
handles.RGBMIPS=rgb;

setappdata(handles.MainFigureSC,'fig',handles.RGBMIPS);
guidata(hObject,handles)
display('image adjusted')

%check_Callback(hObject,handles);


function scalered_Callback(hObject, eventdata, handles)
% hObject    handle to scalered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scalered as scalered3
%        str2double(get(hObject,'String')) returns contents of scalered as a double
handles.proc.scalered = str2num(get(handles.scalered,'String'));

set(handles.autoscale,'Value',0);
handles.justone=1;
guidata(hObject,handles)
check_Callback(hObject, handles);
handles = guidata(hObject);
handles.justone=0;

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function scalered_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scalered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scalegreen_Callback(hObject, eventdata, handles)
% hObject    handle to scalered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scalered as scalered3
%        str2double(get(hObject,'String')) returns contents of scalered as a double
handles.proc.scalegreen = str2num(get(handles.scalegreen,'String'));

set(handles.autoscale,'Value',0);
handles.justone=1;
guidata(hObject,handles)
check_Callback(hObject,handles);
handles.justone=0;

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function scalegreen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scalered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in tiffs.
function tiffs_Callback(hObject, eventdata, handles)
% hObject    handle to tiffs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% scans = str2num(get(handles.startframe,'String'));
% ss = size(handles.inputdata);
% % need to handle normalization better
dataDirectory = handles.readgui.directory;

if(isdir([dataDirectory , '/', 'tiff_stacks']) == 0)
    mkdir ([dataDirectory , '/','tiff_stacks']);
end
scanName = handles.readgui.info.info.scanName;

if (isdir([dataDirectory '/tiff_stacks/', scanName]) == 0)
    mkdir([dataDirectory '/tiff_stacks/', scanName]);
end
%             vals = str2num(get(handles.interpval,'String'));

FilePath = [dataDirectory '/tiff_stacks/', scanName,'/'];

% need to decide whether to use the ones for the current frame (auto or user), or over the
% whole dataset

scr = str2num(get(handles.scalered,'String'));
scg = str2num(get(handles.scalegreen(1,1),'String'));
minr = str2num(get(handles.minr,'String'));
ming = str2num(get(handles.ming,'String'));

handles.writetodesktop=0;
guidata(hObject,handles);
if get(handles.preview_crop,'Value')
    % quick check
    i = max(str2num(get(handles.startframe,'String')));
    
    [red, green] = go_Callback(hObject,handles,i);
    
    if length(size(red))==3
        rgb(:,:,1) = squeeze(scr*(double(max(red(:,:,:),[],2))-minr));
        rgb(:,:,2) = squeeze(scg*(double(max(green(:,:,:),[],2))-ming));
        rgb(:,:,3) = zeros(size(rgb(:,:,1)));
        rgb2(:,:,1) = squeeze(scr*(double(max(red(:,:,:),[],3))-minr));
        rgb2(:,:,2) = squeeze(scg*(double(max(green(:,:,:),[],3))-ming));
        rgb2(:,:,3) = zeros(size(rgb2(:,:,1)));
    else
        rgb(:,:,1) = squeeze(scr*(double(max(red(:,:,:,round(end/2)),[],2))-minr));
        rgb(:,:,2) = squeeze(scg*(double(max(green(:,:,:,round(end/2)),[],2))-ming));
        rgb(:,:,3) = zeros(size(rgb(:,:,1)));
        rgb2(:,:,1) = squeeze(scr*(double(max(red(:,:,:,round(end/2)),[],3))-minr));
        rgb2(:,:,2) = squeeze(scg*(double(max(green(:,:,:,round(end/2)),[],3))-ming));
        rgb2(:,:,3) = zeros(size(rgb2(:,:,1)));
    end
    
    temp11 = squeeze(uint8(256*rgb));
    temp12 =squeeze(uint8(256*rgb2));
    figure(100);
    clf
    imagesc(temp11);
    title('select crop area (lateral - top-left, then bottom-right)');
    [y1, x1] = ginput_col(1);
    hold on
    plot(y1,x1,'.w');
    [y2, x2] = ginput_col(1);
    hold on
    plot(y2,x2,'.w');
    crp0 = floor(y1);
    crp0y = floor(x1);
    crp1 = ceil(y2);
    crp1y = ceil(x2);
    figure(100)
    clf
    imagesc(temp12);
    title('select crop area (depth - left, then right)');
    [z1, y1] = ginput_col(1);
    hold on
    plot(z1,y1,'.w');
    hold on
    [z2, y2] = ginput_col(1);
    plot(z2,y2,'.w');
    crp2 = floor(z1);
    crp3 = ceil(z2);
    clf
    chk='yes';
    clear rgb rgb2
end

disp('Writing tiffs. Please wait.')
separateTiffs = get(handles.SeparateTiffs,'Value');
skewbool = get(handles.Skewcorrection,'Value');
delta = str2num(get(handles.SkewAng,'String'));
m=0;
conversionFactors = [handles.readgui.info.info.cal.ylat; handles.readgui.info.info.cal.zdep; handles.readgui.info.info.cal.xwid];
if length(handles.ss)==3; handles.ss(4) = 1; end

%set cropping margins if not already set
[red, green] = go_Callback(hObject,handles,1);
if ~get(handles.preview_crop,'Value')
    crp0 = 1;
    crp1 = size(red,3);
    crp0y = 1;
    crp1y = size(red,1);
    crp2 = 1;
    crp3 = size(red,2);
end
if crp0==0; crp0 = 1; crp1 = size(red,3); crp0y = 1; crp1y = size(red,1); crp2 = 1; crp3 = size(red,2); end
orient_lat = handles.orient_lat;
if handles.ss(4)==1
    
    fileCounter = 1;
    if ~get(orient_lat, 'Value')
        if separateTiffs
            imgToSave1_HR = [FilePath,'R_unCorrected_' scanName, '.tiff'];
            imgToSave2_HR = [FilePath,'G_unCorrected_' scanName, '.tiff'];
            while (2 == exist(imgToSave1_HR, 'file'))
                imgToSave1_HR = [FilePath,'R_unCorrected_' scanName, '_' num2str(fileCounter) '.tiff'];
                imgToSave2_HR = [FilePath,'G_unCorrected_' scanName, '_' num2str(fileCounter) '.tiff'];
                fileCounter = fileCounter + 1;
            end
        else
            imgToSave_HR = [FilePath,'unCorrected_RGB_' scanName, '.tiff'];
            while (2 == exist(imgToSave_HR, 'file'))
                imgToSave_HR = [FilePath,'unCorrected_RGB_' scanName, '_' num2str(fileCounter) '.tiff'];
                fileCounter = fileCounter + 1;
            end
        end
    else
        if separateTiffs
            imgToSave1_HR = [FilePath,'R_unSkewed_' scanName, '.tiff'];
            imgToSave2_HR = [FilePath,'G_unSkewed_' scanName, '.tiff'];
            while (2 == exist(imgToSave1_HR, 'file'))
                imgToSave1_HR = [FilePath,'R_unSkewed_' scanName, '_' num2str(fileCounter) '.tiff'];
                imgToSave2_HR = [FilePath,'G_unSkewed_' scanName, '_' num2str(fileCounter) '.tiff'];
                fileCounter = fileCounter + 1;
            end
        else
            imgToSave_HR = [FilePath,'unSkewed_RGB_' scanName, '.tiff'];
            while (2 == exist(imgToSave_HR, 'file'))
                imgToSave_HR = [FilePath,'unSkewed_RGB_' scanName, '_' num2str(fileCounter) '.tiff'];
                fileCounter = fileCounter + 1;
            end
        end
    end
    
    red = squeeze(double(red(crp0y:crp1y,crp2:crp3,crp0:crp1)));
    green = squeeze(double(green(crp0y:crp1y,crp2:crp3,crp0:crp1)));
    
    if skewbool
        fileCounter = 1;
        if separateTiffs
            imgToSave1_HR = [FilePath,'R_' scanName, '_R.tiff'];
            imgToSave2_HR = [FilePath,'G_' scanName, '_G.tiff'];
            while (2 == exist(imgToSave1_HR, 'file'))
                imgToSave1_HR = [FilePath,'R_' scanName, '_' num2str(fileCounter) '.tiff'];
                imgToSave2_HR = [FilePath,'G_' scanName, '_' num2str(fileCounter) '.tiff'];
                fileCounter = fileCounter + 1;
            end
        else
            imgToSave_HR = [FilePath,'RGB_' scanName, '.tiff'];
            while (2 == exist(imgToSave_HR, 'file'))
                imgToSave_HR = [FilePath,'RGB_' scanName, '_' num2str(fileCounter) '.tiff'];
                fileCounter = fileCounter + 1;
            end
        end
        % Coordinate System Correction
        red = flip(flip(red,2),3);
        red = permute(red, [3 1 2]);
        RR = imref3d(size(red), conversionFactors(1), conversionFactors(3), conversionFactors(2));
        
        green = flip(flip(green(:,:,:),2),3);
        green = permute(green, [3 1 2]);
        
        RG = imref3d(size(green), conversionFactors(1), conversionFactors(3), conversionFactors(2));
        
        affineMatrix = [1 0 0 0;
            0 1 0 0;
            0 cotd(delta) 1 0;
            0 0 0 1];
        tform = affine3d(affineMatrix);
        
        [red, ~] = imwarp(red, RR, tform);
        [green, ~] = imwarp(green, RG, tform);
    elseif get(handles.orient_lat, 'Value')
        % Coordinate System Correction
        red = flip(flip(red,2),3);
        red = permute(red, [3 1 2]);
        green = flip(flip(green(:,:,:),2),3);
        green = permute(green, [3 1 2]);
    end
    
    
    
    if get(handles.ScaleTiffs,'Value')
        red = uint16(2^16*squeeze(scr*(double(red)-minr)));
        green = uint16(2^16*squeeze(scg*(double(green)-ming)));
    else
        red = uint16(squeeze(double(red)-minr));
        green = uint16(squeeze(double(green)-ming));
    end
    
    
    ss = size(red);
    for i = 1:ss(3)
        if separateTiffs
            imwrite(red(:,:,i), imgToSave1_HR, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
            imwrite(green(:,:,i), imgToSave2_HR, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
        else
            clear rgb
            rgb(:,:,1) = red(:,:,i);
            rgb(:,:,2) = green(:,:,i);
            rgb(:,:,3) = zeros(size(red(:,:,i)));
            imwrite(rgb, imgToSave_HR, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
        end
    end
    
    if get(handles.MIP,'Value')
        tmpr = max(max(red(:,:,:),[],3)-min(red(:,:,:),[],3),[],2);
        [ar, br] = hist(reshape(double(tmpr),[1,prod(size(tmpr))]),200);
        cc1r = br(min(find(ar>0.01*mean(ar))));
        cc2r = br(max(find(ar>0.05*mean(ar))));
        tmpg = max(max(green(:,:,:),[],3)-min(green(:,:,:),[],3),[],2);
        [ag, bg] = hist(reshape(double(tmpg),[1,prod(size(tmpg))]),200);
        cc1g = bg(min(find(ag>0.01*mean(ag))));
        cc2g = bg(max(find(ag>0.05*mean(ag))));
        
        tempr = squeeze(max(red(:,2:end-2,:),[],2))';
        tempr = cat(1,tempr,1000*ones(2,size(red,1)));
        tempr = cat(1,tempr, (squeeze(max(red(:,2:end-2,:),[],3))'));
        tempr = uint8((256/(cc2r-cc1r))*(tempr-cc1r));
        tempg = squeeze(max(green(:,2:end-2,:),[],2))';
        tempg = cat(1,tempg,1000*ones(2,size(green,1)));
        tempg = cat(1,tempg, (squeeze(max(green(:,2:end-2,:),[],3))'));
        tempg = uint8((256/(cc2g-cc1g))*(tempg-cc1g));
        temprg = cat(3,tempr, tempg, zeros(size(tempr)));
        clear tempg tempr tmpr tmpg
        filename =  fullfile(dataDirectory, ['RGB' scanName '_',num2str(handles.readgui.info.info.daq.scanRate),'VPS_',get(handles.readgui.secs,'String'),'secs_MIP.jpg']);
        imwrite(temprg,filename,'jpg')
        clear temprg
    end
else
    clear rgb
    [red, green] = go_Callback(hObject,handles,1);
    tmpr = max(max(red(:,:,:),[],3)-min(red(:,:,:),[],3),[],2);
    [ar, br] = hist(reshape(double(tmpr),[1,prod(size(tmpr))]),200);
    cc1r = br(min(find(ar>0.01*mean(ar))));
    cc2r = br(max(find(ar>0.01*mean(ar))));
    tmpg = max(max(green(:,:,:),[],3)-min(green(:,:,:),[],3),[],2);
    [ag, bg] = hist(reshape(double(tmpg),[1,prod(size(tmpg))]),200);
    cc1g = bg(min(find(ag>0.01*mean(ag))));
    cc2g = bg(max(find(ag>0.01*mean(ag))));
    if get(handles.MIP,'Value')
        disp('Writing MIP movie..')
        filename =  fullfile(dataDirectory, ['RGB' scanName '_',num2str(handles.readgui.info.info.daq.scanRate),'VPS_',get(handles.readgui.secs,'String'),'secs_MIP.avi']);
        vidObj = VideoWriter(filename, 'Uncompressed AVI');
        vidObj.FrameRate = (handles.readgui.info.info.daq.scanRate)/str2num(get(handles.readgui.dsf,'String'));
        open(vidObj);
    end
    clear red green
    
    s4 = handles.ss(4);
    parfor_progress(s4);
   try
        % timepoint = handles.timepoint;
        xg = handles.proc.rgb.xg2;
        xr = handles.proc.rgb.xr2;
        yr = handles.proc.rgb.yr2;
        yg = handles.proc.rgb.yg2;
        zr = handles.proc.rgb.zr2;
        zg = handles.proc.rgb.zg2;
              
        left = floor(min([xg xr]))-1;
        right = floor(size(handles.inputdata,1)-xr);
        
        top = floor(min([zr zg]))-1;
        bot = floor(size(handles.inputdata,2)-max([zr zg]));
        topy = floor(min([yr yg]))-1;
        boty = floor(size(handles.inputdata,3)-max([yr yg]));
        
       % parfor kk =1:s4
        for kk =1:s4
            fileCounter = 1;
            if ~get(orient_lat, 'Value')
                if separateTiffs
                    imgToSave1 = [FilePath,'R_unCorrected_' scanName, '_t', num2str(kk) '.tiff'];
                    imgToSave2 = [FilePath,'G_unCorrected_' scanName, '_t', num2str(kk) '.tiff'];
                    while (2 == exist(imgToSave1, 'file'))
                        imgToSave1 = [FilePath,'R_unCorrected_' scanName, '_' num2str(fileCounter) '_t', num2str(kk) '.tiff'];
                        imgToSave2 = [FilePath,'G_unCorrected_' scanName, '_' num2str(fileCounter) '_t', num2str(kk) '.tiff'];
                        fileCounter = fileCounter + 1;
                    end
                else
                    imgToSave = [FilePath,'unCorrected_RGB_' scanName, '_t', num2str(kk) '.tiff'];
                    while (2 == exist(imgToSave, 'file'))
                        imgToSave = [FilePath,'unCorrected_RGB_' scanName, '_' num2str(fileCounter) '_t', num2str(kk) '.tiff'];
                        fileCounter = fileCounter + 1;
                    end
                end
            else
                if separateTiffs
                    imgToSave1 = [FilePath,'R_unSkewed_' scanName, '_t', num2str(kk) '.tiff'];
                    imgToSave2 = [FilePath,'G_unSkewed_' scanName, '_t', num2str(kk) '.tiff'];
                    while (2 == exist(imgToSave1, 'file'))
                        imgToSave1 = [FilePath,'R_unSkewed_' scanName, '_' num2str(fileCounter) '_t', num2str(kk) '.tiff'];
                        imgToSave2 = [FilePath,'G_unSkewed_' scanName, '_' num2str(fileCounter) '_t', num2str(kk) '.tiff'];
                        fileCounter = fileCounter + 1;
                    end
                else
                    imgToSave = [FilePath,'unSkewed_RGB_' scanName, '_t', num2str(kk) '.tiff'];
                    while (2 == exist(imgToSave, 'file'))
                        imgToSave = [FilePath,'unSkewed_RGB_' scanName, '_' num2str(fileCounter) '_t', num2str(kk) '.tiff'];
                        fileCounter = fileCounter + 1;
                    end
                end
            end
            
            
            % warning - calling Go_Callback within a parfor loop replicates SCAPE_data variable
            % - do not call this function!!
            % [red, green] = go_Callback(hObject,handles,kk);
            
            green = (handles.inputdata(round(xg+[-left:right]),round(zg+[-top:bot]),round(yg+[-topy:boty]),kk));
            red = (handles.inputdata(round(xr+[-left:right]),round(zr+[-top:bot]),round(yr+[-topy:boty]),kk));
            
            %scale gradient because imagesplitter leads to shading in opposite
            %sides of red and green images
            if get(handles.gradient,'Value') == 1
                [grad, ~, ~] = meshgrid(1:size(red,1), 1:size(red,2), 1:size(red,3));
                grad = ((grad-1)./(size(red,1)-1))*2+1;
                grad = flipud(permute(grad, [2 1 3]));
                green = uint16(double(green).*grad);
                red = uint16(double(red).*flipud(grad));
            end
            
            % rotation correction for 1st gen imagesplitter
            if (str2num(get(handles.RotateGreen,'String'))~=0)
                if get(handles.MIPview,'Value')== 2
                    rotval = str2num(get(handles.RotateGreen,'String'));
                    green = imrotate(green,rotval,'bicubic','Crop');
                end
            end
            % scale correction for 1st gen imagesplitter
            if (str2num(get(handles.ScaleG,'String'))~=1)
                scaleval = str2num(get(handles.ScaleG,'String'));
                %             if get(handles.MIPview,'Value')== 2 % scale only xy dimensions
                %                 green = imresize3(green(:,:,iii),[scaleval*size(green,1) size(green,2) size(green,3)]);
                %             else
                %                 green = imresize3(green(:,:,iii),[scaleval*size(green,1) scaleval*size(green,2) size(green,3)]);
                %             end
                green = imresize(green,scaleval, 'bicubic');
                if str2num(get(handles.ScaleG,'String')) <1
                    temp = zeros(size(red));
                    dif1 = size(temp,1) - size(green,1);
                    dif2 =  size(temp,2) - size(green,2);
                    x1 = floor(dif1/2)+1 ;
                    x2 = ceil(dif1/2);
                    y1 = floor(dif2/2)+1;
                    y2 = ceil(dif2/2);
                    temp(x1:end-x2,y1:end-y2,:)=green;
                    green = temp;
                else
                    dif1 = size(green,1) - size(red,1);
                    dif2 =  size(green,2) - size(red,2);
                    x1 = floor(dif1/2)+1; x2 = ceil(dif1/2);
                    y1 = floor(dif2/2)+1; y2 = ceil(dif2/2);
                    green = green(x1:end-x2,y1:end-y2,:);
                end
            end
            
            red = squeeze(double(red(crp0y:crp1y,crp2:crp3,crp0:crp1)));
            green = squeeze(double(green(crp0y:crp1y,crp2:crp3,crp0:crp1)));
            
            
            if skewbool
                fileCounter = 1;
                if separateTiffs
                    imgToSave1 = [FilePath,'R_' scanName, '_t', num2str(kk) '.tiff'];
                    imgToSave2 = [FilePath,'G_' scanName, '_t', num2str(kk) '.tiff'];
                    while (2 == exist(imgToSave1, 'file'))
                        imgToSave1 = [FilePath,'R_' scanName, '_' num2str(fileCounter) '_t', num2str(kk) '.tiff'];
                        imgToSave2 = [FilePath,'G_' scanName, '_' num2str(fileCounter) '_t', num2str(kk) '.tiff'];
                        fileCounter = fileCounter + 1;
                    end
                    
                else
                    imgToSave = [FilePath,'RGB_' scanName, '_t', num2str(kk) '.tiff'];
                    while (2 == exist(imgToSave, 'file'))
                        imgToSave = [FilePath,'RGB' scanName, '_' num2str(fileCounter) '_t', num2str(kk) '.tiff'];
                        fileCounter = fileCounter + 1;
                    end
                end
                % Coordinate System Correction
                red = flip(flip(red(:,:,:),2),3);
                red = permute(red, [3 1 2]);
                RR = imref3d(size(red), conversionFactors(1), conversionFactors(3), conversionFactors(2));
                
                green = flip(flip(green(:,:,:),2),3);
                green = permute(green, [3 1 2]);
                RG = imref3d(size(green), conversionFactors(1), conversionFactors(3), conversionFactors(2));
                
                affineMatrix = [1 0 0 0;
                    0 1 0 0;
                    0 cotd(delta) 1 0;
                    0 0 0 1];
                tform = affine3d(affineMatrix);
                [red , ~] = imwarp(red, RR, tform);
                [green , ~] = imwarp(green, RG, tform);
                
            elseif get(handles.orient_lat, 'Value')
                % Coordinate System Correction
                red = flip(flip(red,2),3);
                red = permute(red, [3 1 2]);
                green = flip(flip(green(:,:,:),2),3);
                green = permute(green, [3 1 2]);
            end
            ss = size(red);
            
            if get(handles.ScaleTiffs,'Value')
                red = uint16(2^16*squeeze(scr*(double(red)-minr)));
                green = uint16(2^16*squeeze(scg*(double(green)-ming)));
            else
                red = uint16(squeeze(double(red)-minr));
                green = uint16(squeeze(double(green)-ming));
            end
            for j = 1:ss(3)
                if separateTiffs
                    imwrite(red(:,:,j), imgToSave1, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                    imwrite(green(:,:,j), imgToSave2, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                else
                    rgb = cat(3, red(:, :, j), green(:, :, j), zeros(size(red(:,:,j))));
                    imwrite(squeeze(rgb), imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                    
                end
            end
            
            if get(handles.MIP,'Value')
                tempr = squeeze(max(red(:,2:end-2,:),[],2))';
                tempr = cat(1,tempr,1000*ones(2,size(red,1)));
                tempr = cat(1,tempr, (squeeze(max(red(:,2:end-2,:),[],3))'));
                tempr = uint8((256/(cc2r-cc1r))*(tempr-cc1r));
                tempg = squeeze(max(green(:,2:end-2,:),[],2))';
                tempg = cat(1,tempg,1000*ones(2,size(green,1)));
                tempg = cat(1,tempg, (squeeze(max(green(:,2:end-2,:),[],3))'));
                tempg = uint8((256/(cc2g-cc1g))*(tempg-cc1g));
                temprg = cat(3,tempr, tempg, zeros(size(tempr)));
                M(kk) = im2frame(temprg);
            end
            
            parfor_progress;
        end
        if get(handles.MIP, 'Value')
            writeVideo(vidObj,M);
            close(vidObj)
        end
        if skewbool
            disp('Skew correction applied')
        end
        
        %toc
    catch
        keyboard
        disp('unsuccessful - if writing MIPs, do not use parfor')
    end
end
disp('Finished writing tiffs');




% --- Executes on slider movement.
function redslide_Callback(hObject, eventdata, handles)
% hObject    handle to redslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

scr = handles.proc.scalered;%
s = 1-get(handles.redslide,'Value');
scr = 1/(s*2/scr);
set(handles.scalered,'String',num2str(scr));
%handles.proc.scalered = scr;
set(handles.autoscale,'Value',0);
handles.justone=1;
guidata(hObject,handles)
check_Callback(hObject, handles);
handles.justone=0;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function redslide_CreateFcn(hObject, ~, ~)
% hObject    handle to redslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function greenslide_Callback(hObject, eventdata, handles)
% hObject    handle to greenslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
scg = handles.proc.scalegreen;%
s = 1-get(handles.greenslide,'Value');
scg = 1/(s*2/scg);
set(handles.scalegreen,'String',num2str(scg));
%handles.proc.scalegreen = scg;
set(handles.autoscale,'Value',0);
handles.justone=1;
guidata(hObject,handles)
check_Callback(hObject,handles);
handles.justone=0;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function greenslide_CreateFcn(hObject, ~, ~)
% hObject    handle to greenslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in autoscale.
function autoscale_Callback(hObject, eventdata, handles)
% hObject    handle to autoscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.autoscale,'Value')
    set(handles.scalegreen,'String',handles.proc.autosc_scalegreen);
    set(handles.greenslide,'Value',0.5);
    handles.proc.scalegreen = handles.proc.autosc_scalegreen;
    set(handles.scalered,'String',handles.proc.autosc_scalered);
    set(handles.redslide,'Value',0.5);
    handles.proc.scalered = handles.proc.autosc_scalered;
    set(handles.minr,'String',handles.proc.autosc_minred);
    set(handles.ming,'String',handles.proc.autosc_minred);
    check_Callback(hObject,handles);
    guidata(hObject,handles)
end
% Hint: get(hObject,'Value') returns toggle state of autoscale



function minr_Callback(hObject, eventdata, handles)
% hObject    handle to minr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minr as scalered3
%        str2double(get(hObject,'String')) returns contents of minr as a double
handles.justone=1;
guidata(hObject,handles)
check_Callback(hObject,handles);
handles.justone=0;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function minr_CreateFcn(hObject, ~, handles)
% hObject    handle to minr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ming_Callback(hObject, eventdata, handles)
% hObject    handle to ming (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ming as scalered3
%        str2double(get(hObject,'String')) returns contents of ming as a double
handles.justone=1;
guidata(hObject,handles)
check_Callback(hObject, handles);
handles.justone=0;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ming_CreateFcn(hObject, ~, ~)
% hObject    handle to ming (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in preview_crop.
function preview_crop_Callback(~, ~, ~)
% hObject    handle to preview_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of preview_crop


% --- Executes on button press in Skewcorrection.
function Skewcorrection_Callback(~, ~, ~)
% hObject    handle to Skewcorrection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Skewcorrection


% --- Executes on button press in MIP.
function MIP_Callback(hObject, eventdata, handles)
% hObject    handle to MIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function SkewAng_Callback(hObject, eventdata, handles)
% hObject    handle to SkewAng (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SkewAng as scalered3
%        str2double(get(hObject,'String')) returns contents of SkewAng as a double


% --- Executes during object creation, after setting all properties.
function SkewAng_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SkewAng (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SeparateTiffs.
function SeparateTiffs_Callback(hObject, eventdata, handles)
% hObject    handle to SeparateTiffs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SeparateTiffs



function RotateGreen_Callback(hObject, eventdata, handles)
% hObject    handle to RotateGreen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RotateGreen as scalered3
%        str2double(get(hObject,'String')) returns contents of RotateGreen as a double
% rgb = getappdata(handles.MainFigureSC,'fig');
% if (str2num(get(handles.RotateGreen,'String'))~=0)
%
%     axis(handles.MainFigureSC);
%     keyboard
%     rotval = str2num(get(handles.RotateGreen,'String'))-handles.rotpastval;
%     if get(handles.MIPview,'Value')== 2
%     rgb(:,:,2) = imrotate(rgb(:,:,2),rotval,'bicubic','Crop');
%     end
%     imagesc(squeeze(uint8((256)*rgb)));
%     axis off
%     handles.rotpastval = str2num(get(handles.RotateGreen,'String'));
%     handles.RGBMIPS = rgb;
%
%     setappdata(handles.MainFigureSC,'fig',rgb);
%     guidata(hObject,handles)
%     display('image adjusted')
%
% end
check_Callback(hObject, handles)
handles.rotpastval = str2num(get(handles.RotateGreen,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function RotateGreen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RotateGreen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PreviewSkew.
function PreviewSkew_Callback(hObject, eventdata, handles)
% hObject    handle to PreviewSkew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.writetodesktop =0;
[red, ~] = go_Callback(hObject,handles,max(str2num(get(handles.startframe,'String'))));
% Coordinate System Correction
conversionFactors = [handles.readgui.info.info.cal.ylat; handles.readgui.info.info.cal.zdep; handles.readgui.info.info.cal.xwid];
red = flip(flip(red,2),3);
red = permute(red, [3 1 2]);
RR = imref3d(size(red), conversionFactors(1), conversionFactors(3), conversionFactors(2));

% Correct for skew
figure(99); subplot(2,1,1)
imagesc([0:size(red, 1)-1]*conversionFactors(3),[0:size(red, 3)-1]*conversionFactors(2), (double(squeeze(max(red(:,:, :), [], 2)))'))
colormap gray; axis image; title('Raw')


delta = str2num(get(handles.SkewAng,'String'));
%delta = pi/2;
affineMatrix = [1 0 0 0;
    0 1 0 0;
    0 cotd(delta) 1 0;
    0 0 0 1];
tform = affine3d(affineMatrix);

[red, ~] = imwarp(red, RR, tform);
figure(99); subplot(2,1,2)
imagesc([0:size(red, 1)-1]*conversionFactors(3),[0:size(red, 3)-1]*conversionFactors(2), (double(squeeze(max(red(:,:, :), [], 2)))'))
colormap gray; axis image;  title('Skew Corrected')
disp('Skew correction applied')


% --- Executes on button press in SaveTransforms.
function SaveTransforms_Callback(hObject, eventdata, handles)
% hObject    handle to SaveTransforms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

transforms.rotateg = str2num(get(handles.RotateGreen,'String'));
transforms.scaleg = str2num(get(handles.ScaleG,'String'));
transforms.scr = str2num(get(handles.scalered,'String'));
transforms.scg = str2num(get(handles.scalegreen(1,1),'String'));
transforms.minr = str2num(get(handles.minr,'String'));
transforms.ming = str2num(get(handles.ming,'String'));
transforms.xg = handles.proc.rgb.xg2;
transforms.xr =handles.proc.rgb.xr2;
transforms.yr = handles.proc.rgb.yr2;
transforms.yg = handles.proc.rgb.yg2;
transforms.zr = handles.proc.rgb.zr2;
transforms.zg = handles.proc.rgb.zg2;

FilePath = fullfile(handles.readgui.directory, ['RGBtransform.mat']);

% checks for repeat file names to prevent overwrite of previous info files
fileCounter = 1;
while (2 == exist(FilePath, 'file'))
    FilePath =  fullfile(handles.readgui.directory, ['RGBtransform_' num2str(fileCounter) '.mat']);
    fileCounter = fileCounter + 1;
end

% save transform info file
save(FilePath, 'transforms')
disp('Transforms saved')


% --- Executes on button press in LoadTransforms.
function LoadTransforms_Callback(hObject, eventdata, handles)
% hObject    handle to LoadTransforms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

presentDirectory = pwd;
cd(handles.readgui.directory);
[fileName, pathName, ~] = uigetfile('*.mat');
cd (pathName)
load(fileName);
cd(presentDirectory)

set(handles.scalered,'String',num2str(transforms.scr));
set(handles.scalegreen,'String',num2str(transforms.scg));
set(handles.minr,'String',num2str(transforms.minr));
set(handles.ming,'String',num2str(transforms.ming));
handles.proc.rgb.xg2 = transforms.xg;
handles.proc.rgb.xr2 = transforms.xr;
handles.proc.rgb.yr2 = transforms.yr;
handles.proc.rgb.yg2 = transforms.yg;
handles.proc.rgb.zr2 = transforms.zr;
handles.proc.rgb.zg2 = transforms.zg;
disp('Transforms applied')


% --- Executes on button press in ScaleTiffs.
function ScaleTiffs_Callback(hObject, eventdata, handles)
% hObject    handle to ScaleTiffs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ScaleTiffs



function ScaleG_Callback(hObject, eventdata, handles)
% hObject    handle to ScaleG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScaleG as scalered3
%        str2double(get(hObject,'String')) returns contents of ScaleG as a double
handles.RGBMIPS = getappdata(handles.MainFigureSC,'fig');

if (str2num(get(handles.ScaleG,'String'))~=0) || (str2num(get(handles.ScaleG,'String'))~=handles.scalepastval)
    
    scaleval = str2num(get(handles.ScaleG,'String'))/handles.scalepastval;
    
    %if (str2num(get(handles.ScaleG,'String'))~=1)
    
    green = handles.RGBMIPS(:,:,2);
    if get(handles.MIPview,'Value')== 2
        green = imresize(green,scaleval);
    else
        green = imresize(green,[(size(green,1)) scaleval*(size(green,2))]);
    end
    
    if str2num(get(handles.ScaleG,'String'))-handles.scalepastval <0
        temp = zeros(size(handles.RGBMIPS(:,:,2)));
        dif1 = size(temp,1) - size(green,1);
        dif2 =  size(temp,2) - size(green,2);
        x1 = floor(dif1/2)+1; x2 = ceil(dif1/2);
        y1 = floor(dif2/2)+1; y2 = ceil(dif2/2);
        
        temp(x1:end-x2,y1:end-y2)=green;
        green = temp;
    else
        dif1 = size(green,1) - size(handles.RGBMIPS(:,:,2),1);
        dif2 =  size(green,2) - size(handles.RGBMIPS(:,:,2),2);
        x1 = floor(dif1/2)+1; x2 = ceil(dif1/2);
        y1 = floor(dif2/2)+1; y2 = ceil(dif2/2);
        green = green(x1:end-x2,y1:end-y2);
    end
    handles.RGBMIPS(:,:,2) = green;
    %end
    
    
    setappdata(handles.MainFigureSC,'fig',handles.RGBMIPS);
    axis(handles.MainFigureSC);
    imagesc(squeeze(uint8((256)*handles.RGBMIPS)));
    axis off
    handles.scalepastval = str2num(get(handles.ScaleG,'String'));
    guidata(hObject,handles)
    display('image adjusted')
end


% --- Executes during object creation, after setting all properties.
function ScaleG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScaleG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scalegreen10_Callback(hObject, eventdata, handles)
% hObject    handle to scalegreen10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scalegreen10 as scalered3
%        str2double(get(hObject,'String')) returns contents of scalegreen10 as a double


% --- Executes during object creation, after setting all properties.
function scalegreen10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scalegreen10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in gradient.
function gradient_Callback(hObject, eventdata, handles)
% hObject    handle to gradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
check_Callback(hObject, handles)


% --- Executes on button press in check2.
function check2_Callback(hObject, eventdata, handles)
% hObject    handle to check2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
check_Callback(hObject, handles)

% --- Executes on button press in go_export.
function go_export_Callback(hObject, eventdata, handles)
% hObject    handle to go_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tic
h = waitbar(0,'Writing to workspace');
handles.writetodesktop=1;
s4 = handles.ss(4);

[redtemp, greentemp] = go_Callback(hObject,handles,1);
red = zeros([size(redtemp) s4]);
green = zeros([size(greentemp) s4]);
red(:,:,:,1) = redtemp;
green(:,:,:,1) = greentemp;

for nn = 2:s4
    [redtemp, greentemp] = go_Callback(hObject,handles,nn);
    red(:,:,:,nn) = redtemp;
    green(:,:,:,nn) = greentemp;
    waitbar(nn/s4,h);
end
assignin('base','red',red);
assignin('base','green',green);
assignin('base','params',handles.proc);
assignin('base','scan_info',handles.readgui);
close(h)
toc

% --- Executes on button press in orient_lat.
function orient_lat_Callback(hObject, eventdata, handles)
% hObject    handle to orient_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of orient_lat
