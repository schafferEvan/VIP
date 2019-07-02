function varargout = SCAPE_read_GUI_Universal8_FLY(varargin)
% SCAPE_read_GUI_Universal8_FLY MATLAB code for SCAPE_read_GUI_Universal8_FLY.fig
%      SCAPE_read_GUI_Universal8_FLY, by itself, creates a new SCAPE_read_GUI_Universal8_FLY or raises the existing_tiff
%      singleton*.
%
%      H = SCAPE_read_GUI_Universal8_FLY returns the handle to a new SCAPE_read_GUI_Universal8_FLY or the handle to
%  java -versio    the existing_tiff singleton*.
%
%      SCAPE_read_GUI_Universal8_FLY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCAPE_read_GUI_Universal8_FLY.M with the given input arguments.
%
%      SCAPE_read_GUI_Universal8_FLY('Property','Value',...) creates a new SCAPE_read_GUI_Universal8_FLY or raises the
%      existing_tiff singleton*.  Starting from the left, property value gridpairs are
%      applied to the GUI before SCAPE_read_GUI_Universal8_FLY_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SCAPE_read_GUI_Universal8_FLY_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SCAPE_read_GUI_Universal8_FLY

% Last Modified by GUIDE v2.5 21-Aug-2018 03:10:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SCAPE_read_GUI_Universal8_FLY_OpeningFcn, ...
    'gui_OutputFcn',  @SCAPE_read_GUI_Universal8_FLY_OutputFcn, ...
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


% --- Executes just before SCAPE_read_GUI_Universal8_FLY is made visible.
function SCAPE_read_GUI_Universal8_FLY_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SCAPE_read_GUI_Universal8_FLY (see VARARGIN)
% set(handles.foldername, 'String', 'F:/');
% handles.initialfolder = '/local_mount/space/enterprise/4/SCAPE_data_backup/';
currentFolder = pwd;
addpath(currentFolder);

handles.initialfolder = 'F:\';
set(handles.foldername, 'String', handles.initialfolder);

% disp('Creating spool file names...')
handles.nooutput = 0;
% Update handles structure
guidata(hObject, handles)


% Creates a list of spool file names that exist for a given run this is important as the camera generates a list that appears unintelligible
clear namesOut tempName
spoolCounter=0;
for i = 1:5000
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
handles.filesToLoad = namesOut;
clear namesOut tempName temp a spoolCounter
guidata(hObject,handles);
% Choose default command line output for SCAPE_read_GUI_Universal8_FLY
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SCAPE_read_GUI_Universal8_FLY wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SCAPE_read_GUI_Universal8_FLY_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
startDirectory = handles.initialfolder;%get(handles.foldername,'String');
directory = [uigetdir(startDirectory,'Select Experiment directory') '/'];
if 0 == directory
    errordlg('You did not select a directory. Please choose again.','yeah')
    return;
end

directory = directory(1:end-1);
set(handles.foldername,'String',directory);

filesInDataDirectory = dir(directory);
nameCounter=0;
for i = 1:length(filesInDataDirectory)
    if filesInDataDirectory(i).isdir ==1 && strcmp(filesInDataDirectory(i).name,'movies')==0
        nameCounter=nameCounter+1;
        runNames{nameCounter} = filesInDataDirectory(i).name;
    end
end

% [SELECTION,OK] = listdlg('ListString',runNames,'ListSize',[400 300]);
handles.directory = directory;
handles.runNames = runNames;
set(handles.listbox1,'Value',1);
set(handles.listbox1,'String',runNames);
[handles] = listbox1_Callback(hObject, eventdata, handles);
guidata(hObject,handles);

% --- Executes on selection change in listbox1.
function [handles] = listbox1_Callback(hObject, ~, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

directory = handles.directory;
filesInDataDirectory = dir(directory);
nameCounter=0;
for i = 1:length(filesInDataDirectory)
    if filesInDataDirectory(i).isdir ==1 && strcmp(filesInDataDirectory(i).name,'movies')==0
        nameCounter=nameCounter+1;
        runNames{nameCounter} = filesInDataDirectory(i).name;
    end
end
set(handles.listbox1,'String',runNames);


qqq = get(handles.listbox1,'Value');
runNames = get(handles.listbox1,'String');
scanName = runNames{qqq};

set(handles.listbox1,'String',runNames);

if strcmp(scanName,'.')==0 && strcmp(scanName,'..')==0
    handles.scanName = scanName;
    infoFilePath = fullfile(handles.directory, [scanName '_info.mat']);
    %     cd (dataDirectory)
    if exist(infoFilePath)
        % load info file
        handles.info = load(infoFilePath);
    end
    % save memory
    if isfield(handles.info.info.daq,'scanWaveform')
        handles.info.info.daq = rmfield(handles.info.info.daq,'scanWaveform');
        handles.info.info.daq.stim = rmfield(handles.info.info.daq.stim,'Vout');
    end
    % Creates file path of zyla metadata file (.ini)
     zylaInfoFilePath = fullfile(handles.directory, scanName, 'acquisitionmetadata.ini');
    if (2 == exist(zylaInfoFilePath))
        FID = fopen(zylaInfoFilePath, 'r');
        zylaMetaData = fread(FID, '*char')';
        fclose(FID);
        % Reads the number of frames in a spool file
        handles.info.zyla.numDepths = str2double(zylaMetaData(12+(strfind(zylaMetaData, 'AOIHeight = ')):(strfind(zylaMetaData, 'AOIWidth = ')-1)));
        handles.info.zyla.numLatPix = round(str2double(zylaMetaData((11+strfind(zylaMetaData, 'AOIStride = ')):(strfind(zylaMetaData, 'PixelEncoding = ')-1)))/2);
        handles.info.zyla.imageBytes = str2double(zylaMetaData((17+strfind(zylaMetaData, 'ImageSizeBytes = ')):(strfind(zylaMetaData, '[multiimage]')-1)));
        handles.info.zyla.numFramesPerSpool = str2double(zylaMetaData((16+strfind(zylaMetaData, 'ImagesPerFile = ')):end));
        startIndex = strfind(zylaMetaData, 'ImageSizeBytes') + length('ImageSizeBytes = ');
        handles.info.zyla.ImageSize = str2double(zylaMetaData(startIndex : startIndex + 7));
    else
        disp('No metadata file found for this run');
    end
    
    %generate string of info:
    clear infotext
%     keyboard
    if (handles.info.info.daq.scanRate ~= handles.info.info.camera.framerate/handles.info.info.daq.pixelsPerLine)
        handles.info.info.daq.scanRate = handles.info.info.camera.framerate/handles.info.info.daq.pixelsPerLine;
        info = handles.info.info;
        save(infoFilePath, '-v7.3', 'info')
    end
    
    infotext{1} = sprintf('Scan Rate: %.2f VPS',handles.info.info.daq.scanRate);
    infotext{2} = sprintf('Scan Duration: %.2f secs',handles.info.info.daq.scanLength);
    infotext{3} = sprintf('X-steps: %d, X-range: %.4f V',handles.info.info.daq.pixelsPerLine,handles.info.info.daq.scanAngle);
    infotext{4} = sprintf('Camera Frame Rate: %.2f Hz',handles.info.info.daq.pixelrate);
    infotext{5} = sprintf('Camera height: %d, Width: %d',handles.info.zyla.numDepths,handles.info.zyla.numLatPix);
    
    if(isfield(handles.info.info,'sawtooth'))
        infotext{6} = sprintf('Sawtooth: Yes');
        set(handles.saw,'Value',1);
    else
        if datenum(handles.info.info.scanStartTimeApprox) > 7.3621e+05
            infotext{6} = sprintf('Sawtooth: Probably yes');
            set(handles.saw,'Value',1);
        else
            infotext{6} = sprintf('Sawtooth: No');
            set(handles.saw,'Value',0);
        end
    end
    
    % calibration factors
%         scandate = datenum(handles.info.info.scanStartTimeApprox);
%         calpath = 'C:\Users\Goat\Documents\MATLAB\SCAPE\Load Code\';
%         if exist([calpath 'SCAPE_calibrations.txt'],'file')==2
%             fid = fopen([calpath 'SCAPE_calibrations.txt'],'r');
%             fgetl(fid);
%             %         fgetl(fid);
%             m=1;
%             a{m} = fgetl(fid);
%             while a{m}~=-1
%                 m=m+1;
%                 a{m}= fgetl(fid);
%                 p = strfind(a{m},'	');
%                 try sdate = datenum(a{m}(1:p(1)-1));
%                     edate = datenum(a{m}(p(1)+1:p(2)-1));
%                     if scandate>=sdate && scandate<=edate
%                         handles.info.info.cal = [];
%                         handles.info.info.cal.ylat = str2num(a{m}(p(2)+1:p(3)-1));
%                         handles.info.info.cal.zdep = str2num(a{m}(p(3)+1:p(4)-1));
%                         handles.info.info.cal.xfact = str2num(a{m}(p(4)+1:end));
%                         handles.info.info.cal.xwid = handles.info.info.cal.xfact*handles.info.info.daq.scanAngle/handles.info.info.daq.pixelsPerLine;
%                         a{m} = -1;
%                         success = 1;
%                     end
%                 catch
%                     success = 0;
%                 end
%             end
%             fclose(fid)
%             if success ==1
%                 infotext{8} = sprintf('Calibration factors: y: %.2f  x: %.1f z: %.1f um/pix',handles.info.info.cal.ylat, handles.info.info.cal.xwid,handles.info.info.cal.zdep);
%                 infotext{9} = sprintf('Field of view: y: %.2f  x: %.1f z: %.1f um',handles.info.info.cal.ylat*handles.info.zyla.numLatPix, handles.info.info.cal.xwid*handles.info.info.daq.pixelsPerLine,handles.info.zyla.numDepths*handles.info.info.cal.zdep);
%             else
%                 disp('! calibration factors for this DATE not found');
%                 handles.info.info.cal.ylat = 1.23 ;%1.24;
%                 handles.info.info.cal.zdep = 	1.03; %1.22;
%                 handles.info.info.cal.xfact = 	300;% 327.28;
%                 handles.info.info.cal.xwid = handles.info.info.cal.xfact*handles.info.info.daq.scanAngle/handles.info.info.daq.pixelsPerLine;
%                 handles.info.info.cal2 = sprintf('Get recent vals! Using xwid=%0.2f ylat=1.24 zdep=1.22', handles.info.info.cal.xwid);
%                 infotext{8} = sprintf('Calibration: %s',handles.info.info.cal2);
%             end
%         else
%             disp('! calibration factor FILE not found');
%     
%             handles.info.info.cal.ylat = 1.23 ;%1.24;
%             handles.info.info.cal.zdep = 	1.03; %1.22;
%             handles.info.info.cal.xfact = 	300;% 327.28;
%             handles.info.info.cal.xwid = handles.info.info.cal.xfact*handles.info.info.daq.scanAngle/handles.info.info.daq.pixelsPerLine;
%             handles.info.info.cal2 = sprintf('Get recent vals! Using xwid=%0.2f ylat=1.24 zdep=1.22', handles.info.info.cal.xwid);
%             infotext{8} = sprintf('Calibration: %s',handles.info.info.cal2);
%         end
%     
    
    try
        handles.info.info.cal.ylat = handles.info.info.GUIcalFactors.y_umPerPix;
        handles.info.info.cal.zdep = handles.info.info.GUIcalFactors.z_umPerPix;
        handles.info.info.cal.xfact = handles.info.info.GUIcalFactors.xK_umPerVolt;
        handles.info.info.cal.xwid = handles.info.info.cal.xfact*handles.info.info.daq.scanAngle/handles.info.info.daq.pixelsPerLine;
        infotext{8} = sprintf('Calibration factors: y: %.2f  x: %.1f z: %.1f um/pix',handles.info.info.cal.ylat, handles.info.info.cal.xwid,handles.info.info.cal.zdep);
        infotext{9} = sprintf('Field of view: y: %.2f  x: %.1f z: %.1f um',handles.info.info.cal.ylat*handles.info.zyla.numLatPix, handles.info.info.cal.xwid*handles.info.info.daq.pixelsPerLine,handles.info.zyla.numDepths*handles.info.info.cal.zdep);
        
    catch
        disp('! calibration factors not found, using 1.23 y um/pix, 1.03 z um/pix, and 300 x um/volt');
        handles.info.info.cal.ylat = 1.37;
        handles.info.info.cal.zdep = 1.15;
        handles.info.info.cal.xfact = 317.9;
        handles.info.info.cal.xwid = handles.info.info.cal.xfact*handles.info.info.daq.scanAngle/handles.info.info.daq.pixelsPerLine;
        infotext{8} = sprintf('Calibration factors: y: %.2f  x: %.1f z: %.1f um/pix',handles.info.info.cal.ylat, handles.info.info.cal.xwid,handles.info.info.cal.zdep);
        infotext{9} = sprintf('Field of view: y: %.2f  x: %.1f z: %.1f um',handles.info.info.cal.ylat*handles.info.zyla.numLatPix, handles.info.info.cal.xwid*handles.info.info.daq.pixelsPerLine,handles.info.zyla.numDepths*handles.info.info.cal.zdep);
        
    end
    
    infotext{7} = sprintf('Scan date and time: %s',handles.info.info.scanStartTimeApprox);
    infotext{10} = sprintf('Est whole dataset file: %.1f Mb',handles.info.info.camera.xROI*handles.info.info.camera.yROI*str2num(handles.info.info.camera.kineticSeriesLength)/1000000);
    infotext{11} = sprintf('(camfile) Est whole dataset file: %.1f Mb',handles.info.zyla.imageBytes*str2num(handles.info.info.camera.kineticSeriesLength)/1000000);
    infotext{12} = sprintf('Experiment Notes: %s',mat2str(handles.info.info.experiment_notes(1,:)));
    for kk = 2:size(handles.info.info.experiment_notes,1)
        infotext{kk+11} = sprintf('%s',mat2str(handles.info.info.experiment_notes(kk,:)));
    end
    numSpoolfilesPerVolume = handles.info.info.daq.pixelsPerLine/handles.info.zyla.numFramesPerSpool;
    numActual_spoolfiles = length(dir([handles.directory,'/',handles.scanName,'/*.dat']));
    numspools = numSpoolfilesPerVolume*handles.info.info.daq.scanRate*handles.info.info.daq.scanLength;
    
    if numspools>numActual_spoolfiles
        infotext{11} = sprintf('POSSIBLE FAILED RUN (%g files out of %g)',numActual_spoolfiles, ceil(numspools));
        handles.pf = 'pf';
    else
        handles.pf = '';
    end
    
    m=1;
    existing{1} = 'no TIFF stacks found';
    %  infotext{9} = 'existing_tiff outputs found:';
    if isdir([handles.directory,'/tiff_stacks'])
        tiffs = dir([handles.directory,'/tiff_stacks']);
        for i = 1:length(tiffs)
            if strfind(tiffs(i).name,scanName)
                
                numtiffs = dir([handles.directory,'/tiff_stacks/',tiffs(i).name,'/*.tiff']);
                existing{m} = sprintf('%s  (%d files)',tiffs(i).name, length(numtiffs));
                m=m+1;
                
            end
        end
    end
    set(handles.existing_tiff,'String',existing);
    
    m=1;
    set(handles.existing_avi,'Value',1);
    existing_avi{1} = 'none';
    avis = dir([handles.directory,'/*.avi']);
    pics = dir([handles.directory,'/*.jpg']);
    for i = 1:length(avis)
        if strfind(avis(i).name,scanName)
            existing_avi{m} = sprintf('%s',avis(i).name);
            m=m+1;
        end
    end
    for i = 1:length(pics)
        if strfind(pics(i).name,scanName)
            existing_avi{m} = sprintf('%s',pics(i).name);
            m=m+1;
        end
    end
    set(handles.existing_avi,'String',existing_avi);
    set(handles.text2,'String',infotext);
    
    directory = handles.directory;
    filesInDataDirectory = dir(directory);
    nameCounter=0;
    for i = 1:length(filesInDataDirectory)
        if filesInDataDirectory(i).isdir ==1 && strcmp(filesInDataDirectory(i).name,'movies')==0
            nameCounter=nameCounter+1;
            runNames{nameCounter} = filesInDataDirectory(i).name;
        end
    end
    
    % [SELECTION,OK] = listdlg('ListString',runNames,'ListSize',[400 300]);
    temp = get(handles.listbox1,'Value');
    handles.runNames = runNames;
    handles.scanName = scanName;
    set(handles.listbox1,'Value',1);
    set(handles.listbox1,'String',runNames);
    set(handles.listbox1,'Value',temp);
    guidata(hObject,handles);
end

% --- Executes on button press in preview.
function [framey] = preview_Callback(hObject, eventdata, handles)
% hObject    handle to preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes during object creation, after setting all properties.
[handles] = listbox1_Callback(hObject, eventdata, handles);
numSpoolfilesPerVolume = ceil(handles.info.info.daq.pixelsPerLine/handles.info.zyla.numFramesPerSpool);
numFramesInThisSpool = handles.info.zyla.numFramesPerSpool;
carrydata = [];
frameCounterStart = 1;
h = waitbar(0,'hang on');
depths = handles.info.zyla.numDepths+2;
lateral = handles.info.zyla.numLatPix;

for spoolFileCounter = 1:numSpoolfilesPerVolume+1
    waitbar(spoolFileCounter/numSpoolfilesPerVolume,h);
    
    loadfile = handles.filesToLoad{spoolFileCounter};
    filePath = fullfile(handles.directory, handles.scanName, loadfile);
    FID = fopen(filePath, 'r');
    % Read all data within spool file as a 1D stream of integers and add it to carrydata
    rawData = [carrydata; fread(FID, 'uint16=>uint16')];
    %     numPixelsToReshape = numRows * numColumns * numFramesInThisSpool;
    % close file
    fclose(FID);
    frameCounter  = numFramesInThisSpool*(spoolFileCounter);
    frameCounterStart = numFramesInThisSpool*(spoolFileCounter-1) + 1;
    numPixelsToReshape=(depths*lateral*numFramesInThisSpool);
    dataTemp(:,:,frameCounterStart:frameCounter) = reshape(rawData(1:numPixelsToReshape), [lateral,depths, numFramesInThisSpool]);
    disp('.')
    carrydata = rawData((numPixelsToReshape+1):end);
    clear rawData;
    if frameCounter>handles.info.info.daq.pixelsPerLine
        SCAPE_data(:,:,1:handles.info.info.daq.pixelsPerLine) = dataTemp(:,:,1:handles.info.info.daq.pixelsPerLine);
    end
end
close(h)
figure(100)
subplot(2,1,1)
imagesc(squeeze(max(SCAPE_data(1:end-10,1:end-2,:),[],2))');
axis image
xlabel('Y (pixels)')
ylabel('X (pixels)')
title('X-Y MIP')
subplot(2,1,2)
imagesc((squeeze(max(SCAPE_data(1:end-10,1:end-2,:),[],3))'));
axis image
colormap gray
xlabel('Y (pixels)')
ylabel('Z (pixels)')
title('Y-Z MIP')
framey =  SCAPE_data;

% --- Executes on selection change in existing_avi.
function existing_avi_Callback(~, ~, handles)
% hObject    handle to existing_avi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns existing_avi contents as cell array
%        contents{get(hObject,'Value')} returns selected item from existing_avi
vids = get(handles.existing_avi,'String');
vidchose = get(handles.existing_avi,'Value');
if strcmp(char(vids(vidchose)),'none')==0
    [~,~,ext] = fileparts(char(vids(vidchose)));
    if ext == '.jpg'
        h = figure(101);
        clf
        image1 = imread([handles.directory,'/',char(vids(vidchose))]);
        imagesc(image1)
        pos = get(gca,'Position');
        set(gca,'Position',[0.01 pos(2) 0.99 pos(4)]);
        titley = strrep(char(vids(vidchose)),'_','-');
        axis tight off image
        title(titley)
        colormap gray
    else
        vidObj = VideoReader([handles.directory,'/',char(vids(vidchose))]);
        % hh = warndlg('STOP?','stop')
        h = figure(101);
        clf
        vidFrame = readFrame(vidObj);
        imagesc(vidFrame);
        pos = get(gca,'Position');
        set(gca,'Position',[0.01 pos(2) 0.99 pos(4)]);
        % currAxes = axes;
        titley = strrep(char(vids(vidchose)),'_','-');
        while hasFrame(vidObj) && ishandle(h)
            vidFrame = readFrame(vidObj);
            if ishandle(h)==0
                break
            else
                imagesc(vidFrame);
                axis tight off image
                title(titley)
                pause(1/vidObj.FrameRate)
            end
        end
    end
end



function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saw.
function saw_Callback(hObject, eventdata, handles)
% hObject    handle to saw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saw



function dsf_Callback(hObject, eventdata, handles)
% hObject    handle to dsf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dsf as text
%        str2double(get(hObject,'String')) returns contents of dsf as a double


% --- Executes during object creation, after setting all properties.
function dsf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dsf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function foldername_Callback(hObject, eventdata, handles)
% hObject    handle to foldername (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of foldername as text
%        str2double(get(hObject,'String')) returns contents of foldername as a double


% --- Executes during object creation, after setting all properties.
function foldername_CreateFcn(hObject, eventdata, handles)
% hObject    handle to foldername (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on button press in incdiff.
function incdiff_Callback(hObject, eventdata, handles)
% hObject    handle to incdiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of incdiff


% --- Executes on button press in tiffstack.
function tiffstack_Callback(hObject, eventdata, handles)
% hObject    handle to tiffstack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% msgbox('doesn''t work yet, sorry - export to workspace and use tiff_write_GUI');
handles.nooutput = 1;
guidata(hObject,handles);
if strcmp(get(handles.secs,'String'),'all')==0
    %     [button] = questdlg('Total duration is not set to ALL - continue?','Yes','Woops - ok, stop');
    disp('WARNING: Total duration is not set to ALL ');
end
button = 'Yes';
if strcmp(button,'Yes')
    [SELECTION,OK] = listdlg('ListString',handles.runNames,'ListSize',[400 300],'InitialValue',get(handles.listbox1,'Value'));
    
    for i = 1:length(SELECTION)
        if strcmp(handles.runNames(SELECTION(i)),'.')==0 && strcmp(handles.runNames(SELECTION(i)),'..')==0
            set(handles.listbox1,'Value',SELECTION(i));
            [handles] = listbox1_Callback(hObject, eventdata, handles);
            handles.tiff = 1;
            handles.nooutput = 1;
            guidata(hObject,handles);
            loaddata_Callback(hObject, eventdata, handles)
            [handles] = listbox1_Callback(hObject, eventdata, handles)
        end
    end
    
    
    handles.nooutput = 0;
    handles.tiff = 0;
    guidata(hObject,handles);
end
handles.nooutput = 0;
handles.tiff = 0;
guidata(hObject,handles);

% --- Executes on button press in makeavi.
function makeavi_Callback(hObject, eventdata, handles)
% hObject    handle to makeavi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nooutput = 1;
guidata(hObject,handles);
if strcmp(get(handles.secs,'String'),'all')==0
    %     [button] = questdlg('Total duration is not set to ALL - continue?','Yes','Woops - ok, stop');
    disp('WARNING: Total duration is not set to ALL ');
end
button = 'Yes';
handles.tiff = 0;
set(handles.makeMIPMovie_checkbox, 'Value', 1);
if strcmp(button,'Yes')
    [SELECTION,OK] = listdlg('ListString',handles.runNames,'ListSize',[400 300],'InitialValue',get(handles.listbox1,'Value'));
    for i = 1:length(SELECTION)
        if strcmp(handles.runNames(SELECTION(i)),'.')==0 && strcmp(handles.runNames(SELECTION(i)),'..')==0
            set(handles.listbox1,'Value',SELECTION(i));
            [handles] = listbox1_Callback(hObject, eventdata, handles)
            handles.nooutput = 1;
            guidata(hObject,handles);
            filename =  fullfile(handles.directory, [handles.scanName '__dsf',num2str(get(handles.dsf,'String')), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',handles.pf,'_quikmovie.avi']);
            if exist(filename)
                disp(sprintf('Movie already exists - skipping (%s)',handles.scanName));
            else
                loaddata_Callback(hObject, eventdata, handles)
                [handles] = listbox1_Callback(hObject, eventdata, handles)
            end
        end
    end
    handles.nooutput = 0;
    guidata(hObject,handles);
end
handles.nooutput = 0;
guidata(hObject,handles);




% --- Executes on button press in playmovie.
function playmovie_Callback(hObject, eventdata, handles)
% hObject    handle to playmovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of playmovie


% --- Executes on button press in playAVI.
function playAVI_Callback(hObject, eventdata, handles)
% hObject    handle to playAVI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of playAVI

function secs_Callback(hObject, eventdata, handles)
% hObject    handle to secs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of secs as text
%        str2double(get(hObject,'String')) returns contents of secs as a double


% --- Executes during object creation, after setting all properties.
function secs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to secs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in existing_tiff.
function existing_tiff_Callback(hObject, eventdata, handles)
% hObject    handle to existing_tiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns existing_tiff contents as cell array
%        contents{get(hObject,'Value')} returns selected item from existing_tiff


% --- Executes during object creation, after setting all properties.
function existing_tiff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to existing_tiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in loaddata.
function loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to loaddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hh2 = waitbar(1,'Working');

[handles] = listbox1_Callback(hObject, eventdata, handles);
dsf = str2num(get(handles.dsf,'String'));
saw = get(handles.saw,'Value');


numSecondsToLoad = get(handles.secs,'String');
if strcmp(numSecondsToLoad,'all')
    numSecondsToLoad = handles.info.info.daq.scanLength;
else
    numSecondsToLoad = str2num(numSecondsToLoad);
end

% Scale data set size by spatial bin factor
scanRate = handles.info.info.daq.scanRate;               % Volumetric Scan Rate (VPS)
numFrames = handles.info.info.daq.pixelrate*handles.info.info.daq.scanLength;% length(info.daq.scanWaveform);  % Total number of frames acquired
numScans = handles.info.info.daq.numberOfScans;          % Total number of volumes acquired
framesPerScan = handles.info.info.daq.pixelsPerLine;%floor(info.camera.framerate/info.daq.scanRate);        % Number of frames per volume

%kpedit - if stage scanning instead of glavo scanning, framesPerScan=2, so
%reset number of volumes to be 1 and for all frames to be read in
%sequence
if framesPerScan == 2
    framesPerScan = round(numFrames);
    numScans = round(numFrames);
    scanRate = 1/handles.info.info.daq.scanLength;
end


numFramesPerSpool = handles.info.zyla.numFramesPerSpool;
numDepths = handles.info.zyla.numDepths;
numLatPix = handles.info.zyla.numLatPix;


% numSpoolfilesPerVolume = ceil(handles.info.info.daq.pixelsPerLine/handles.info.zyla.numFramesPerSpool);
numSpoolfilesPerVolume = ceil(framesPerScan/numFramesPerSpool);
numSpoolfilesPerVolume1 = (framesPerScan/numFramesPerSpool);

numActual_spoolfiles = length(dir([handles.directory,'/',handles.scanName,'/*.dat']));
if numSpoolfilesPerVolume1*scanRate*numSecondsToLoad>numActual_spoolfiles
    numSecondsToLoad = floor(numActual_spoolfiles/(numSpoolfilesPerVolume1*scanRate))-1;
    disp(sprintf('%s: Not enough spool files - possible failed run, loading max available (%d secs)', handles.scanName, numSecondsToLoad));
    pf = 'pf';
else
    pf = '';
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
        tempName = tempName(find(tempName ~= ' '));
        tempName = [tempName 'spool.dat'];
        namesOut{spoolCounter} = tempName;
    end
    namesOut = [{'0000000000spool.dat'} namesOut];
    handles.filesToLoad = namesOut;
    clear namesOut tempName temp a spoolCounter
    guidata(hObject,handles);
end



% add extra columns for 2 buffer rows
if (2 == handles.info.info.camera.binFactor)
    numColumns = numDepths + 1;
    numRows = handles.info.zyla.ImageSize / 2 / numColumns;
    if(mod(numRows, 1) ~= 0)
        numColumns = numDepths+2;
        numRows = handles.info.zyla.ImageSize/2/numColumns;
    end
else
    % KPEDIT_ TEMP 3/22/18 -FIX THIS
    numColumns = numDepths + 2;
    numRows = numLatPix;
end

close(hh2)

if get(handles.cropload,'Value') == 1
    [framey] = preview_Callback(hObject, eventdata, handles);
    figure(102)
    subplot(2,1,1)
    imagesc(log(double(squeeze(max(framey(1:end-10,1:end-2,:),[],2))')));
    axis image
    xlabel('Y (pixels)')
    ylabel('X (pixels)')
    title('X-Y MIP (LOG)')
    subplot(2,1,2)
    imagesc(log(double(squeeze(max(framey(1:end-10,1:end-2,:),[],3))')));
    colormap gray
    xlabel('Y (pixels)')
    ylabel('Z (pixels)')
    title('Y-Z MIP')
    colormap jet
    %     imagesc(squeeze(max(framey,[],2)
    title('Crop in Z+Y: click upper-left and bottom right of ROI')
    clear framey;
    [YY xx] = ginput(2);
    Ycrop(1) = ceil(min(YY));
    Ycrop(2) = floor(max(YY));
    Xcrop(1) = ceil(min(xx));
    Xcrop(2) = floor(max(xx));
    %     cropping = 1
else
    %     cropping = 0;
    Ycrop = [1 numRows-10];
    Xcrop = [1 numColumns-3];
    xx = [1 handles.info.info.daq.pixelsPerLine];
end
handles.info.info.Ycrop = Ycrop;
handles.info.info.Xcrop = Xcrop;

% pre-loaded spool files (add check that there are enough)
filesToLoad = handles.filesToLoad;

% Figures out index for first frame in each volume
frame_starts = 1:framesPerScan:framesPerScan*numScans;

% Figures out which volumes to load
if (length(numSecondsToLoad) == 1)
    volumesToLoad = 1:dsf:round(numSecondsToLoad*scanRate);
else
    volumesToLoad = round(numSecondsToLoad(1)*scanRate):dsf:round(numSecondsToLoad(2)*scanRate);
end

if isempty(volumesToLoad)==0
    
    % Figures out first frame in each volume that we PLAN TO LOAD
    
    startFrameInVolumeToLoad = frame_starts(volumesToLoad);
    
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
    
    if get(handles.makeDouble_checkbox, 'Value')
        if volumesToLoad == 1
            SCAPE_data = zeros(Ycrop(2)-Ycrop(1)+1, Xcrop(2)-Xcrop(1)+1, numFramesPerSpool,  numFileToLoad, 'double');
        else
            SCAPE_data = zeros(Ycrop(2)-Ycrop(1)+1, Xcrop(2)-Xcrop(1)+1, numFramesPerSpool,  numFileToLoad, 'double');
        end
    else
        
        if volumesToLoad == 1
            SCAPE_data = zeros(Ycrop(2)-Ycrop(1)+1, Xcrop(2)-Xcrop(1)+1, numFramesPerSpool,  numFileToLoad, 'uint16');
        else
            SCAPE_data = zeros(Ycrop(2)-Ycrop(1)+1, Xcrop(2)-Xcrop(1)+1, numFramesPerSpool,  numFileToLoad, 'uint16');
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
            filePath = fullfile(handles.directory, handles.scanName, loadfile);
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
            if get(handles.makeDouble_checkbox,'Value')
                SCAPE_data(:,:,:,spoolFileCounter) = double(rawData(Ycrop(1):Ycrop(2),Xcrop(1):Xcrop(2),:));
            else
                SCAPE_data(:,:,:,spoolFileCounter) = rawData(Ycrop(1):Ycrop(2),Xcrop(1):Xcrop(2),:);
            end
            
        end
    else
        
        parfor_progress(numFileToLoad);
        parfor spoolFileCounter = 1:numFileToLoad
            loadfile = filesToLoad{SpoolToLoad(spoolFileCounter)};
            filePath = fullfile(handles.directory, handles.scanName, loadfile);
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
            if get(handles.makeDouble_checkbox,'Value')
                SCAPE_data(:,:,:,spoolFileCounter) = double(rawData(Ycrop(1):Ycrop(2),Xcrop(1):Xcrop(2),:));
            else
                SCAPE_data(:,:,:,spoolFileCounter) = rawData(Ycrop(1):Ycrop(2),Xcrop(1):Xcrop(2),:);
            end
            parfor_progress;
        end
        parfor_progress(0);
    end
    
    if volumesToLoad == 1
        SCAPE_data = reshape(SCAPE_data,[[Ycrop(2)-Ycrop(1)]+1, Xcrop(2)-Xcrop(1)+1, numFramesPerSpool*numFileToLoad]);
        SCAPE_data = reshape(SCAPE_data(:,:,index),[[Ycrop(2)-Ycrop(1)]+1, Xcrop(2)-Xcrop(1)+1, framesPerScan-2, numVolumeToLoad]);
    else
        SCAPE_data = reshape(SCAPE_data,[[Ycrop(2)-Ycrop(1)]+1, Xcrop(2)-Xcrop(1)+1, numFramesPerSpool*numFileToLoad]);
        SCAPE_data = reshape(SCAPE_data(:,:,index),[[Ycrop(2)-Ycrop(1)]+1, Xcrop(2)-Xcrop(1)+1, framesPerScan-2, numVolumeToLoad]);
    end
    
    toc
    % load waitbar
    if strcmp(handles.pf,'pf') && strcmp(get(handles.secs,'String'),'all');
        yy = [sprintf('(pf: %d secs)',numSecondsToLoad)];
    else
        yy = '';
    end
    
    if handles.nooutput == 0;
        assignin('base','SCAPE_data',SCAPE_data);
        handles.info.info.listbox = get(handles.text2,'String');
        assignin('base','info',handles.info.info);
        assignin('base','zyla',handles.info.zyla);
        h = waitbar(1,'SCAPE data written to workspace');
        pause(2);
        close(h)
    else
        if get(handles.makeMIPMovie_checkbox, 'Value') 
        h = waitbar(1,'Writing movie');
        end
    end
    % Flips even frames to account for bidirectional scan. Also waits until second volume is loaded and flipped to introduce the pixel shift.
    if (saw == 0 && (round(dsf/2)== dsf/2)==0)
        SCAPE_data(:, :, :, 1:2:end) = flipdim(SCAPE_data(:, :, :, 1:2:end), 3);
    end
    
    % close(h)
    clear M
    ss = size(SCAPE_data);
    % Create Preview movie/jpg
    if get(handles.makeMIPMovie_checkbox, 'Value') 
        if get(handles.makeMIPMovie_checkbox, 'Value') && (length(ss)==3 || ss(4) == 1)
            filename =  fullfile(handles.directory, [handles.scanName '__HIRES.jpg']);
            if(isdir([handles.directory,'/tiff_stacks']) == 0)
                mkdir([handles.directory,'/tiff_stacks']);
            end
            imgToSave_topMIP = fullfile(handles.directory,'/tiff_stacks/', [handles.scanName '_topMIP_',pf,'_quickMIP.tiff']);
            imgToSave_sideMIP = fullfile(handles.directory,'/tiff_stacks/', [handles.scanName '_sideMIP_', pf,'_quickMIP.tiff']);
                
            tempy = squeeze(max(SCAPE_data(:,3:end-3,:),[],2))';
            imwrite(tempy, imgToSave_topMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');             
            tempy = cat(1,tempy,1000*ones(2,ss(1)));
            tempy = cat(1,tempy, (squeeze(max(SCAPE_data(:,3:end-3,:),[],3))'));
            tmp = max(SCAPE_data(:,3:end-3,:),[],3);

            imwrite(tmp', imgToSave_sideMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');

            [a, b] = hist(reshape(double(tmp),[1,prod(size(tmp))]),200);
            cc1 = b(min(find(a>0.01*mean(a))));
            cc2 = b(max(find(a>0.01*mean(a))));
            tempy = uint8((2^8/(cc2-cc1))*(tempy-cc1));
            imwrite(tempy,filename,'jpg')
            
        else
            if get(handles.incdiff,'Value')==1
                diffy =2; stt = 7;
            else
                diffy = 1;
                stt = 1;
            end
            
            for pp = 1:diffy
                if pp==1
                    filename =  fullfile(handles.directory, [handles.scanName '__dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quikmovie.avi']);
                    imgToSave_topMIP = fullfile(handles.directory,'/tiff_stacks/', [handles.scanName '_topMIP_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP.tiff']);
                    imgToSave_sideMIP = fullfile(handles.directory,'/tiff_stacks/', [handles.scanName '_sideMIP_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP.tiff']);
                
                else
                    filename =  fullfile(handles.directory, [handles.scanName '__DIFF_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'quikmovie.avi']);
                    imgToSave_topMIP = fullfile(handles.directory,'/tiff_stacks/', [handles.scanName '_topMIP_DIFF_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP.tiff']);
                    imgToSave_sideMIP = fullfile(handles.directory,'/tiff_stacks/', [handles.scanName '_sideMIP_DIFF_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP.tiff']);
                
                end
                if exist(filename) 
                    disp('Movie already exists - skipping');
                elseif ~get(handles.makeMIPMovie_checkbox, 'Value') 
                    disp('Movie step skipped!');
                else
                    clear tempy M
                    if pp==1
                        tmp = max(max(SCAPE_data(:,:,:,1),[],4),[],3);
                    else
                        tmp = max(max(max(SCAPE_data(:,:,:,1:end/5),[],4)-min(SCAPE_data(:,:,:,1:end/5),[],4),[],4),[],3);
                    end
                    [a b] = hist(reshape(double(tmp),[1,prod(size(tmp))]),200);
                    cc1 = b(min(find(a>0.01*mean(a))))*.95;
                    cc2 = b(max(find(a>0.01*mean(a))))*.95;
                    if get(handles.makeMIPMovie_checkbox, 'Value')
                        disp('Writing MIP movie..')
                        vidObj = VideoWriter(filename, 'Uncompressed AVI');
                        vidObj.FrameRate = handles.info.info.daq.scanRate/dsf;
                        open(vidObj);
                        map = gray(256);
                    end
                    for volumeCounter = stt:ss(4)
                        numst = 2; % stitch every nth volume
                        clear tempy
                        if pp==1
                            tempy = squeeze(max(SCAPE_data(:,10:end-10,:,volumeCounter),[],2))';
                            imwrite(uint16(tempy), imgToSave_topMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                            
                            tempy = cat(1,tempy,1000*ones(2,ss(1)));
                           
                            tempy = cat(1,tempy, (squeeze(max(SCAPE_data(:,10:end-10,:,volumeCounter),[],3))'));
                            imwrite(uint16(squeeze(max(SCAPE_data(:,3:end-3,:,volumeCounter),[],3))'), imgToSave_sideMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');

                            % disp(volumeCounter)
                            
                        else
                            tempy = squeeze(max(SCAPE_data(:,10:end-10,:,volumeCounter)-SCAPE_data(:,10:end-10,:,volumeCounter-6),[],2))';
                            imwrite(uint16(tempy), imgToSave_topMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                            
                            tempy = cat(1,tempy,1000*ones(2,ss(1)));
                            tempy = cat(1,tempy, (squeeze(max(SCAPE_data(:,10:end-10,:,volumeCounter)-SCAPE_data(:,10:end-10,:,volumeCounter-6),[],3))'));
                            imwrite(uint16(squeeze(max(SCAPE_data(:,3:end-3,:,volumeCounter),[],3))'), imgToSave_sideMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');                        
                        
                        end
                        if get(handles.makeMIPMovie_checkbox, 'Value')
                            M(volumeCounter-stt+1) = im2frame(uint8((256/(cc2-cc1))*(tempy-cc1)),map);
                        end
                    end
                    if get(handles.makeMIPMovie_checkbox, 'Value')
                        writeVideo(vidObj,M);
                        close(vidObj)
                        close(h)
                    end
               
                    disp('done')
                end
            end
        end
    end
    
    if isfield(handles,'tiff')
        if handles.tiff == 1
            if(isdir([handles.directory,'/tiff_stacks']) == 0)
                mkdir([handles.directory,'/tiff_stacks']);
            end
            
            if(isdir([handles.directory,'/tiff_stacks/',handles.scanName]) == 0)
                mkdir([handles.directory,'/tiff_stacks/',handles.scanName]);
            end
            
            filePath = [handles.directory,'/tiff_stacks/',handles.scanName,'/'];
            %filePath = [handles.directory,'/tiff_stacks/'];

            tic
            disp('Start Saving Scan');
            ss = size(SCAPE_data);
            handles.ss = ss;
            conversionFactors = [handles.info.info.cal.ylat; handles.info.info.cal.zdep; handles.info.info.cal.xwid];
            if numel(ss)== 3
                if ~get(handles.orient_lat, 'Value')
                imgToSave = [filePath, 'unCorrected_' handles.scanName '.tiff'];
                fileCounter = 1;
                while (2 == exist(imgToSave, 'file'))
                    imgToSave = [filePath, 'unCorrected_' handles.scanName '_' num2str(fileCounter) '.tiff'];
                    fileCounter = fileCounter + 1;
                end
                else
                imgToSave = [filePath, 'unSkewed_' handles.scanName '.tiff'];
                fileCounter = 1;
                while (2 == exist(imgToSave, 'file'))
                    imgToSave = [filePath, 'unSkewed_' handles.scanName '_' num2str(fileCounter) '.tiff'];
                    fileCounter = fileCounter + 1;
                end
                end
              
                
                if get(handles.checkbox2,'Value')
                    
                    % Coordinate System Correction
                    SCAPE_data = flip(flip(SCAPE_data(:,:,:),2),3);
                    SCAPE_data = permute(SCAPE_data, [3 1 2]);

                    R = imref3d(size(SCAPE_data), conversionFactors(1), conversionFactors(3), conversionFactors(2));
                    
                    %Correct for skew
                    delta = str2num(get(handles.SkewAng,'String'));
                    
                    affineMatrix = [1 0 0 0;
                        0 1 0 0;
                        0 cotd(delta) 1 0;
                        0 0 0 1];
                    tform = affine3d(affineMatrix);
                    [SCAPE_data, ~] = imwarp(SCAPE_data, R, tform);
                    imgToSave = [filePath, handles.scanName '.tiff'];
                    disp('Skew correction applied')
                elseif get(handles.orient_lat, 'Value')
                     % Coordinate System Correction
                    SCAPE_data = flip(flip(SCAPE_data(:,:,:),2),3);
                    SCAPE_data = permute(SCAPE_data, [3 1 2]);
                end

                for j = 1:size(SCAPE_data, 3)
                    imwrite(uint16(SCAPE_data(:, 2:end-1, j))-105, imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                end
            else
                parfor_progress(ss(4));
                
                parfor i  = 1:ss(4)

                    if ~get(handles.orient_lat, 'Value')
                    imgToSave = [filePath, 'unCorrected_' handles.scanName '_t' num2str(i) '.tiff'];
                    fileCounter = 1;
                    while (2 == exist(imgToSave, 'file'))
                        imgToSave = [filePath, 'unCorrected_' handles.scanName '_' num2str(fileCounter) '_t' num2str(i) '.tiff'];
                        fileCounter = fileCounter + 1;
                    end
                    else
                    imgToSave = [filePath, 'unSkewed_' handles.scanName '_t' num2str(i) '.tiff'];
                    fileCounter = 1;
                    while (2 == exist(imgToSave, 'file'))
                        imgToSave = [filePath, 'unSkewed_' handles.scanName '_' num2str(fileCounter) '_t' num2str(i) '.tiff'];
                        fileCounter = fileCounter + 1;
                    end 
                    end
                    
                    temp = SCAPE_data(:,:,:,i);
                    if get(handles.checkbox2,'Value')
                        % Coordinate System Correction
                        temp = flip(flip(temp,2),3);
                        temp = permute(temp, [3 1 2]);
                        R = imref3d(size(temp), conversionFactors(1), conversionFactors(3), conversionFactors(2));
                        
                        % Correct for skew
                        %                     figure(1)
                        %                     imagesc([0:size(SCAPE_data, 1)-1]*conversionFactors(3), ...
                        %                         [0:size(SCAPE_data, 3)-1]*conversionFactors(2), squeeze(max(SCAPE_data(:, 300:350, :), [], 2))')
                        %                     axis image
                        delta = str2num(get(handles.SkewAng,'String'));
                        %delta = pi/2;
                        affineMatrix = [1 0 0 0;
                            0 1 0 0;
                            0 cotd(delta) 1 0;
                            0 0 0 1];
                        tform = affine3d(affineMatrix);
                        [temp, ~] = imwarp(temp, R, tform);
                        % temp = log(double(temp)+1)*2^16/log(2^16)
                        imgToSave = [filePath, handles.scanName '_' num2str(i) '.tiff'];
                    elseif get(handles.orient_lat, 'Value')
                        % Coordinate System Correction
                        temp = flip(flip(temp(:,:,:),2),3);
                        temp = permute(temp, [3 1 2]);
                    end
                    for j = 1:size(temp, 3)
                        imwrite(uint16(temp(:, 1:end, j))-105, imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                    end
                    parfor_progress;
                end
                parfor_progress(0);
                
            end
            disp('Finished Saving Scan');
            toc
            if get(handles.OpeninFiji,'Value')
                addpath('/local_mount/space/enterprise/4/fiji-linux64-20170530/Fiji.app/scripts/')
                javaaddpath '/usr/local/MATLAB/R2016a/java/mij.jar'
                javaaddpath '/usr/local/MATLAB/R2016a/java/ij.jar'
                Miji;
                
                if numel(ss) == 3
                    numtiffs = 1;
                    if get(handles.checkbox2, 'Value')
                        imgToSaveA = [filePath, handles.scanName '.tiff'];
                    else
                        imgToSaveA = [filePath, 'unCorrected_' handles.scanName '.tiff'];
                    end
                else
                    numtiffs = ss(4);
                    if get(handles.checkbox2, 'Value')
                        imgToSaveA = [filePath, handles.scanName '_t1.tiff'];
                    else
                        imgToSaveA = [filePath, 'unCorrected_' handles.scanName '_t1.tiff'];
                    end
                end
                
                info = imfinfo(imgToSaveA);
                z = length(info);
                
                % command1 = strcat('open=[',directory,'] number=', num2str(numtiffs), ' file=.tiff sort');
                % MIJ.run('Image Sequence...', command1)
                
                command1 = strcat('open=[',imgToSaveA,'] number=', num2str(numtiffs), ' starting=1 increment=1 scale=100 file=.tiff');
                MIJ.run('Hypervolume Opener', command1)
                
                command2 = strcat('order=xyczt(default) channels=1 slices=', num2str(z), ' frames=', num2str(numtiffs), ' display=Color');
                MIJ.run('Stack to Hyperstack...', command2);
                
            end
            
        end
    end
    
    [handles] = listbox1_Callback(hObject, eventdata, handles);
    handles.nooutput = 0; 
    guidata(hObject,handles);
    clear SCAPE_data;
    
else
    disp(sprintf('NOTHING TO LOAD: %s',handles.info.info.scanName));
end


% --- Executes during object creation, after setting all properties.
function existing_avi_CreateFcn(hObject, ~, ~)
% hObject    handle to existing_avi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cropload.
function cropload_Callback(hObject, eventdata, handles)
% hObject    handle to cropload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cropload


% --- Executes on button press in checkbox2.
function checkbox2_Callback(~, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% SkewCorrection
% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in RGBmerge.
function RGBmerge_Callback(~, ~, handles)
% hObject    handle to RGBmerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Call Splitcols

% SCAPE_splitcolsV6B(handles)
% if exist('SCAPE_data')
% SCAPE_splitcolsV7(evalin('base', 'SCAPE_data'),evalin('base', 'info'),handles.directory)
% else

SCAPE_splitcolsV9b(handles)

% --- Executes on button press in makeMIPMovie_checkbox.
function makeMIPMovie_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to makeMIPMovie_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of makeMIPMovie_checkbox


% --- Executes on button press in makeDouble_checkbox.
function makeDouble_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to makeDouble_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of makeDouble_checkbox



function SkewAng_Callback(hObject, eventdata, handles)
% hObject    handle to SkewAng (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SkewAng as text
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


% --- Executes on button press in PreviewSkew.
function PreviewSkew_Callback(hObject, eventdata, handles)
% hObject    handle to PreviewSkew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles] = listbox1_Callback(hObject, eventdata, handles);
numSpoolfilesPerVolume = ceil(handles.info.info.daq.pixelsPerLine/handles.info.zyla.numFramesPerSpool);
numFramesInThisSpool = handles.info.zyla.numFramesPerSpool;
carrydata = [];
frameCounterStart = 1;
h = waitbar(0,'hang on');
depths = handles.info.zyla.numDepths+2;
lateral = handles.info.zyla.numLatPix;
spoolFileCounter = 1;
loadfile = handles.filesToLoad{spoolFileCounter};
filePath = fullfile(handles.directory, handles.scanName, loadfile);
FID = fopen(filePath, 'r');
while (spoolFileCounter <= 2*numSpoolfilesPerVolume+2 && FID~=-1)
    
    waitbar(spoolFileCounter/numSpoolfilesPerVolume,h);
    
    % Read all data within spool file as a 1D stream of integers and add it to carrydata
    rawData = [carrydata; fread(FID, 'uint16=>uint16')];
    %     numPixelsToReshape = numRows * numColumns * numFramesInThisSpool;
    % close file
    fclose(FID);
    frameCounter  = numFramesInThisSpool*(spoolFileCounter);
    frameCounterStart = numFramesInThisSpool*(spoolFileCounter-1) + 1;
    numPixelsToReshape=(depths*lateral*numFramesInThisSpool);
    dataTemp(:,:,frameCounterStart:frameCounter) = reshape(rawData(1:numPixelsToReshape), [lateral,depths, numFramesInThisSpool]);
    carrydata = rawData((numPixelsToReshape+1):end);
    clear rawData;
    if frameCounter>handles.info.info.daq.pixelsPerLine
        SCAPE_data(:,:,1:handles.info.info.daq.pixelsPerLine) = dataTemp(:,:,1:handles.info.info.daq.pixelsPerLine);
    end
    
    % for next load
    
    loadfile = handles.filesToLoad{spoolFileCounter+1};
    filePath = fullfile(handles.directory, handles.scanName, loadfile);
    FID = fopen(filePath, 'r');
    spoolFileCounter = spoolFileCounter + 1;
end
close(h)


%Skew and dsiplay
conversionFactors = [handles.info.info.cal.ylat; handles.info.info.cal.zdep; handles.info.info.cal.xwid];

SCAPE_data = flip(flip(SCAPE_data(:,:,:),2),3);
SCAPE_data = permute(SCAPE_data, [3 1 2]);
R = imref3d(size(SCAPE_data), conversionFactors(1), conversionFactors(3), conversionFactors(2));

figure(99); subplot(2,2,1)
imagesc([0:size(SCAPE_data, 1)-1]*conversionFactors(3),[0:size(SCAPE_data, 3)-1]*conversionFactors(2), log(double(squeeze(max(SCAPE_data(:,round(end/3):round(end/2), :), [], 2))))')
colormap gray; axis image; title('Raw')
subplot(2,2,2)
imagesc([0:size(SCAPE_data, 2)-1]*conversionFactors(1),[0:size(SCAPE_data, 1)-1]*conversionFactors(3), log(double(squeeze(max(SCAPE_data(:,:,round(end/3):round(2*end/3)), [], 3)))))
colormap gray; axis image;


delta = str2num(get(handles.SkewAng,'String'));

affineMatrix = [1 0 0 0;
    0 1 0 0;
    0 cotd(delta) 1 0;
    0 0 0 1];
tform = affine3d(affineMatrix);
[SCAPE_data, ~] = imwarp(SCAPE_data, R, tform);

figure(99); subplot(2,2,3)
imagesc([0:size(SCAPE_data, 1)-1]*conversionFactors(3),[0:size(SCAPE_data, 3)-1]*conversionFactors(2), log(double(squeeze(max(SCAPE_data(:,round(end/3):round(end/2), :), [], 2))))')
colormap gray; axis image; title('Skew Corrected')
subplot(2,2,4)
imagesc([0:size(SCAPE_data, 2)-1]*conversionFactors(1),[0:size(SCAPE_data, 1)-1]*conversionFactors(3), log(double(squeeze(max(SCAPE_data(:,:,round(end/3):round(2*end/3)), [], 3)))))
colormap gray; axis image;

disp('Skew correction applied')


% --- Executes on button press in OpeninFiji.
function OpeninFiji_Callback(hObject, eventdata, handles)
% hObject    handle to OpeninFiji (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OpeninFiji


% --- Executes on button press in GoFiji.
function GoFiji_Callback(hObject, eventdata, handles)
% hObject    handle to GoFiji (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

addpath('/local_mount/space/enterprise/4/fiji-linux64-20170530/Fiji.app/scripts/')
javaaddpath '/usr/local/MATLAB/R2016a/java/mij.jar'
javaaddpath '/usr/local/MATLAB/R2016a/java/ij.jar'
Miji;

qqq = get(handles.existing_tiff,'Value');
tiffNames = get(handles.existing_tiff,'String');
tiffName = tiffNames{qqq};

ind = strfind(tiffName, ' ');
tiffName = tiffName(1:ind-1);

dir1 = dir([handles.directory,'/tiff_stacks/',tiffName, '/*tiff']);
numtiffs = length(dir1);

directory = fullfile(handles.directory, 'tiff_stacks', tiffName, dir1(1).name);
info = imfinfo([directory]);
z = length(info);

command1 = strcat('open=[',directory,'] number=', num2str(numtiffs), ' file=.tiff sort');
MIJ.run('Image Sequence...', command1)


% command1 = strcat('open=[',directory,'] number=', num2str(numtiffs), ' starting=1 increment=1 scale=100 file=.tiff');
% 
% MIJ.run('Hypervolume Opener', command1)

command2 = strcat('order=xyczt(default) channels=1 slices=', num2str(z), ' frames=', num2str(numtiffs), ' display=Color');

MIJ.run('Stack to Hyperstack...', command2);





% --- Executes on button press in quickviewer.
function quickviewer_Callback(hObject, eventdata, handles)
% hObject    handle to quickviewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
loaddata_Callback(hObject, eventdata, handles)
SCAPE_quickview_VV_KP4(evalin('base', 'SCAPE_data'),evalin('base', 'info'),0)


% --- Executes on button press in orient_lat.
function orient_lat_Callback(hObject, eventdata, handles)
% hObject    handle to orient_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of orient_lat


function blocksize_Callback(hObject, eventdata, handles)
% hObject    handle to blocksize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blocksize as text
%        str2double(get(hObject,'String')) returns contents of blocksize as a double


% --- Executes during object creation, after setting all properties.
function blocksize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blocksize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dualtiffstack.
function dualtiffstack_Callback(hObject, eventdata, handles)
% hObject    handle to dualtiffstack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nooutput = 1;
guidata(hObject,handles);
if strcmp(get(handles.secs,'String'),'all')==0
    %     [button] = questdlg('Total duration is not set to ALL - continue?','Yes','Woops - ok, stop');
    disp('WARNING: Total duration is not set to ALL ');
end
button = 'Yes';
transFilePath = fullfile(handles.directory, 'RGBtransform.mat');

    load(transFilePath)
    
    handles.proc.rgb.xg2 = transforms.xg;
    handles.proc.rgb.xr2 = transforms.xr;
    handles.proc.rgb.yr2 = transforms.yr;
    handles.proc.rgb.yg2 = transforms.yg;
    handles.proc.rgb.zr2 = transforms.zr;
    handles.proc.rgb.zg2 = transforms.zg;
    
    
if strcmp(button,'Yes')
    [SELECTION,OK] = listdlg('ListString',handles.runNames,'ListSize',[400 300],'InitialValue',get(handles.listbox1,'Value'));
    
    for i = 1:length(SELECTION)
        if strcmp(handles.runNames(SELECTION(i)),'.')==0 && strcmp(handles.runNames(SELECTION(i)),'..')==0
            set(handles.listbox1,'Value',SELECTION(i));
            [handles] = listbox1_Callback(hObject, eventdata, handles);
            handles.tiff = 0;
            handles.nooutput = 1;
            guidata(hObject,handles);
            loaddatablock_Callback(hObject, eventdata, handles)
            [handles] = listbox1_Callback(hObject, eventdata, handles)
        end
    end
    
    
    handles.nooutput = 0;
    handles.tiff = 0;
    guidata(hObject,handles);
end
handles.nooutput = 0;
handles.tiff = 0;
guidata(hObject,handles);

function loaddatablock_Callback(hObject, eventdata, handles)
hh2 = waitbar(1,'Working');

[handles] = listbox1_Callback(hObject, eventdata, handles);
dsf = str2num(get(handles.dsf,'String'));
saw = get(handles.saw,'Value');

xg = handles.proc.rgb.xg2;
xr = handles.proc.rgb.xr2;
yr = handles.proc.rgb.yr2;
yg = handles.proc.rgb.yg2;
zr = handles.proc.rgb.zr2;
zg = handles.proc.rgb.zg2;

numSecondsToLoad = get(handles.secs,'String');
if strcmp(numSecondsToLoad,'all')
    numSecondsToLoad = handles.info.info.daq.scanLength;
else
    numSecondsToLoad = str2num(numSecondsToLoad);
end

% Scale data set size by spatial bin factor
scanRate = handles.info.info.daq.scanRate;               % Volumetric Scan Rate (VPS)
numFrames = handles.info.info.daq.pixelrate*handles.info.info.daq.scanLength;% length(info.daq.scanWaveform);  % Total number of frames acquired
numScans = handles.info.info.daq.numberOfScans;          % Total number of volumes acquired
framesPerScan = handles.info.info.daq.pixelsPerLine;%floor(info.camera.framerate/info.daq.scanRate);        % Number of frames per volume

%kpedit - if stage scanning instead of glavo scanning, framesPerScan=2, so
%reset number of volumes to be 1 and for all frames to be read in
%sequence
if framesPerScan == 2
    framesPerScan = round(numFrames);
    numScans = round(numFrames);
    scanRate = 1/handles.info.info.daq.scanLength;
end


numFramesPerSpool = handles.info.zyla.numFramesPerSpool;
numDepths = handles.info.zyla.numDepths;
numLatPix = handles.info.zyla.numLatPix;


% numSpoolfilesPerVolume = ceil(handles.info.info.daq.pixelsPerLine/handles.info.zyla.numFramesPerSpool);
numSpoolfilesPerVolume = ceil(framesPerScan/numFramesPerSpool);
numSpoolfilesPerVolume1 = (framesPerScan/numFramesPerSpool);

numActual_spoolfiles = length(dir([handles.directory,'/',handles.scanName,'/*.dat']));
if numSpoolfilesPerVolume1*scanRate*numSecondsToLoad>numActual_spoolfiles
    numSecondsToLoad = floor(numActual_spoolfiles/(numSpoolfilesPerVolume1*scanRate))-1;
    disp(sprintf('%s: Not enough spool files - possible failed run, loading max available (%d secs)', handles.scanName, numSecondsToLoad));
    pf = 'pf';
else
    pf = '';
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
        tempName = tempName(find(tempName ~= ' '));
        tempName = [tempName 'spool.dat'];
        namesOut{spoolCounter} = tempName;
    end
    namesOut = [{'0000000000spool.dat'} namesOut];
    handles.filesToLoad = namesOut;
    clear namesOut tempName temp a spoolCounter
    guidata(hObject,handles);
end



% add extra columns for 2 buffer rows
if (2 == handles.info.info.camera.binFactor)
    numColumns = numDepths + 1;
    numRows = handles.info.zyla.ImageSize / 2 / numColumns;
    if(mod(numRows, 1) ~= 0)
        numColumns = numDepths+2;
        numRows = handles.info.zyla.ImageSize/2/numColumns;
    end
else
    % KPEDIT_ TEMP 3/22/18 -FIX THIS
    numColumns = numDepths + 2;
    numRows = numLatPix;
end

close(hh2)

% if get(handles.cropload,'Value') == 1
%     [framey] = preview_Callback(hObject, eventdata, handles);
%     figure(102)
%     subplot(2,1,1)
%     imagesc(log(double(squeeze(max(framey(1:end-10,1:end-2,:),[],2))')));
%     axis image
%     xlabel('Y (pixels)')
%     ylabel('X (pixels)')
%     title('X-Y MIP (LOG)')
%     subplot(2,1,2)
%     imagesc(log(double(squeeze(max(framey(1:end-10,1:end-2,:),[],3))')));
%     colormap gray
%     xlabel('Y (pixels)')
%     ylabel('Z (pixels)')
%     title('Y-Z MIP')
%     colormap jet
%     %     imagesc(squeeze(max(framey,[],2)
%     title('Crop in Z+Y: click upper-left and bottom right of ROI')
%     clear framey;
%     [YY xx] = ginput(2);
%     Ycrop(1) = ceil(min(YY));
%     Ycrop(2) = floor(max(YY));
%     Xcrop(1) = ceil(min(xx));
%     Xcrop(2) = floor(max(xx));
%     %     cropping = 1
% else
%     %     cropping = 0;
%     Ycrop = [1 numRows-10];
%     Xcrop = [1 numColumns-3];
%     xx = [1 handles.info.info.daq.pixelsPerLine];
% end
% handles.info.info.Ycrop = Ycrop;
% handles.info.info.Xcrop = Xcrop;

left = floor(min([xg xr]))-1;
right = floor(numRows-10-xr);

top = floor(min([zr zg]))-1;
bot = floor(numColumns-3-max([zr zg]));

Ycrop_g = round(xg+[-left right]);
Xcrop_g = round(zg+[-top bot]);

Ycrop_r = round(xr+[-left right]);
Xcrop_r = round(zr+[-top bot]);

handles.info.info.Ycrop_g = Ycrop_g;
handles.info.info.Xcrop_g = Xcrop_g;

handles.info.info.Ycrop_r = Ycrop_r;
handles.info.info.Xcrop_r = Xcrop_r;

Ycrop = [1 numRows-10];
Xcrop = [1 numColumns-3];

% pre-loaded spool files (add check that there are enough)
filesToLoad = handles.filesToLoad;

% Figures out index for first frame in each volume
frame_starts = 1:framesPerScan:framesPerScan*numScans;

% Figures out which volumes to load
if (length(numSecondsToLoad) == 1)
    volumesToLoad = 1:dsf:round(numSecondsToLoad*scanRate);
else
    volumesToLoad = round(numSecondsToLoad(1)*scanRate):dsf:round(numSecondsToLoad(2)*scanRate);
end

if isempty(volumesToLoad)==0
    
    % Figures out first frame in each volume that we PLAN TO LOAD
    
    startFrameInVolumeToLoad = frame_starts(volumesToLoad);
    
    disp('Pre-allocating memory, please wait ...')
    
    % Creating the index of every single frames.
    numVolumeToLoad = length(startFrameInVolumeToLoad);
    endFrameInVolumeToLoad = startFrameInVolumeToLoad+framesPerScan-3;
    index_vol = zeros(framesPerScan-2,numVolumeToLoad);
    SpoolToLoad_vol = zeros(numSpoolfilesPerVolume+1,numVolumeToLoad);
    for i = 1:numVolumeToLoad
        index_vol(:,i) = startFrameInVolumeToLoad(i):1:endFrameInVolumeToLoad(i);
        SpoolToLoad_vol(1:length(unique(ceil(index_vol(:,i) /numFramesPerSpool))),i) = unique(ceil(index_vol(:,i) /numFramesPerSpool));

    end
    index = reshape(index_vol,[1,((framesPerScan-2)*numVolumeToLoad)]);
    
    %         SpoolToLoad = unique([ceil(startFrameInVolumeToLoad/numFramesPerSpool),ceil(endFrameInVolumeToLoad/numFramesPerSpool)]);
    SpoolToLoad = unique(ceil(index/numFramesPerSpool));
    numFileToLoad = length(SpoolToLoad);
    numFrameToLoad = numVolumeToLoad*framesPerScan;
    
    numPixelsToReshape = numRows * numColumns * numFramesPerSpool;
    
    startSpoolInVolumeToLoad = ceil(startFrameInVolumeToLoad/numFramesPerSpool);
    startFrameInSpoolToLoad = mod(startFrameInVolumeToLoad,numFramesPerSpool);
    
    % Load all data into this variable (lateral, depth, scan, time)
     SCAPE_data_g = zeros(Ycrop_g(2)-Ycrop_g(1)+1, Xcrop_g(2)-Xcrop_g(1)+1, numFramesPerSpool,  numSpoolfilesPerVolume+1, 'uint16');
     SCAPE_data_r = zeros(Ycrop_r(2)-Ycrop_r(1)+1, Xcrop_r(2)-Xcrop_r(1)+1,  numFramesPerSpool,  numSpoolfilesPerVolume+1, 'uint16');
     SCAPE_data = zeros(Ycrop(2)-Ycrop(1)+1, Xcrop(2)-Xcrop(1)+1, numFramesPerSpool,  numSpoolfilesPerVolume+1, 'double');
     SCAPE_data_temp = zeros(Ycrop_r(2)-Ycrop_r(1)+1, Xcrop_r(2)-Xcrop_r(1)+1,framesPerScan-2);

     SCAPE_data_temp = flip(flip(SCAPE_data_temp,2),3);
     SCAPE_data_temp = permute(SCAPE_data_temp, [3 1 2]);
    
    frameLoadInSpool = zeros(numFramesPerSpool,numFileToLoad);
    for i = 1:numFileToLoad
        frameLoadInSpool(:,i) = ((SpoolToLoad(i)-1)*numFramesPerSpool+1):((SpoolToLoad(i)-1)*numFramesPerSpool+numFramesPerSpool);
    end
    frameLoadInSpool = reshape(frameLoadInSpool,[1,numFramesPerSpool*numFileToLoad]);
    index = find(ismember(frameLoadInSpool,index));
    
    for i = 1:numVolumeToLoad
        index_vol(:,i) = find(ismember(frameLoadInSpool,index_vol(:,i)));
    end
    
    disp('done with memory allocation')
    % calculate the skew-correction
    conversionFactors = [handles.info.info.cal.ylat; handles.info.info.cal.zdep; handles.info.info.cal.xwid];
    R = imref3d(size(SCAPE_data_temp), conversionFactors(1), conversionFactors(3), conversionFactors(2));
    
    %Correct for skew
    delta = str2num(get(handles.SkewAng,'String'));
    
    affineMatrix = [1 0 0 0;
        0 1 0 0;
        0 cotd(delta) 1 0;
        0 0 0 1];
    tform = affine3d(affineMatrix);
    
    if(isdir([handles.directory,'/tiff_stacks']) == 0)
        mkdir([handles.directory,'/tiff_stacks']);
    end
    
    if(isdir([handles.directory,'/tiff_stacks/',handles.scanName]) == 0)
        mkdir([handles.directory,'/tiff_stacks/',handles.scanName]);
    end
    
    filePath_tiff = [handles.directory,'/tiff_stacks/',handles.scanName,'/'];
    
    tic
            
        parfor_progress(numFileToLoad);
    parfor volumeCounter = 1:length(volumesToLoad)
             SCAPE_data_temp = zeros(Ycrop(2)-Ycrop(1)+1, Xcrop(2)-Xcrop(1)+1, numFramesPerSpool,  numSpoolfilesPerVolume+1, 'double');

%         volumeCounter
%         for spoolFileCounter = 1:size(SpoolToLoad_vol(:,volumeCounter),1)
        for spoolFileCounter = 1:length(nonzeros(SpoolToLoad_vol(:,volumeCounter)))
            loadfile = filesToLoad{SpoolToLoad_vol(spoolFileCounter,volumeCounter)};
            filePath = fullfile(handles.directory, handles.scanName, loadfile);
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
            if get(handles.makeDouble_checkbox,'Value')
                SCAPE_data_temp(:,:,:,spoolFileCounter) = double(rawData(Ycrop(1):Ycrop(2),Xcrop(1):Xcrop(2),:));
            else
                SCAPE_data_temp(:,:,:,spoolFileCounter) = rawData(Ycrop(1):Ycrop(2),Xcrop(1):Xcrop(2),:);
            end
%             
        end
        parfor_progress;

        SCAPE_data_rh = reshape(SCAPE_data_temp,[[Ycrop(2)-Ycrop(1)]+1, Xcrop(2)-Xcrop(1)+1, numFramesPerSpool*(numSpoolfilesPerVolume+1)]);
        SCAPE_data_rh = reshape(SCAPE_data_rh(:,:,index_vol(:,volumeCounter)-(SpoolToLoad_vol(1,volumeCounter)-1)*numFramesPerSpool),[[Ycrop(2)-Ycrop(1)]+1, Xcrop(2)-Xcrop(1)+1, framesPerScan-2]);
        SCAPE_data_g = SCAPE_data_rh(Ycrop_g(1):Ycrop_g(2),Xcrop_g(1):Xcrop_g(2),:);
        SCAPE_data_r = SCAPE_data_rh(Ycrop_r(1):Ycrop_r(2),Xcrop_r(1):Xcrop_r(2),:);
        SCAPE_data_g = flip(flip(SCAPE_data_g,2),3);
        SCAPE_data_g = permute(SCAPE_data_g, [3 1 2]);
        SCAPE_data_r = flip(flip(SCAPE_data_r,2),3);
        SCAPE_data_r = permute(SCAPE_data_r, [3 1 2]);
%         R = imref3d(size(SCAPE_data_r), conversionFactors(1), conversionFactors(3), conversionFactors(2));
%         
%         delta = str2num(get(handles.SkewAng,'String'));
%         %delta = pi/2;
%         affineMatrix = [1 0 0 0;
%             0 1 0 0;
%             0 cotd(delta) 1 0;
%             0 0 0 1];
%         tform = affine3d(affineMatrix);
%         [temp, ~] = imwarp(temp, R, tform);
        
        [SCAPE_data_g, ~] = imwarp(SCAPE_data_g, R, tform);
        [SCAPE_data_r, ~] = imwarp(SCAPE_data_r, R, tform);

        imgToSave_g = [filePath_tiff, 'G_' handles.scanName '_' num2str(volumesToLoad(volumeCounter)) '.tiff'];
        imgToSave_r = [filePath_tiff, 'R_' handles.scanName '_' num2str(volumesToLoad(volumeCounter)) '.tiff'];
        for j = 1:size(SCAPE_data_g, 3)
            imwrite(uint16(SCAPE_data_g(:, :, j))-105, imgToSave_g, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
            imwrite(uint16(SCAPE_data_r(:, :, j))-105, imgToSave_r, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
        end
    end
    parfor_progress(0);
    toc
    % load waitbar
    if strcmp(handles.pf,'pf') && strcmp(get(handles.secs,'String'),'all');
        yy = [sprintf('(pf: %d secs)',numSecondsToLoad)];
    else
        yy = '';
    end
    
    if handles.nooutput == 0;
        assignin('base','SCAPE_data',SCAPE_data);
        handles.info.info.listbox = get(handles.text2,'String');
        assignin('base','info',handles.info.info);
        assignin('base','zyla',handles.info.zyla);
        h = waitbar(1,'SCAPE data written to workspace');
        pause(2);
        close(h)
    else
        if get(handles.makeMIPMovie_checkbox, 'Value') 
        h = waitbar(1,'Writing movie');
        end
    end
    % Flips even frames to account for bidirectional scan. Also waits until second volume is loaded and flipped to introduce the pixel shift.
    if (saw == 0 && (round(dsf/2)== dsf/2)==0)
        SCAPE_data(:, :, :, 1:2:end) = flipdim(SCAPE_data(:, :, :, 1:2:end), 3);
    end
    
    % close(h)
    clear M
    ss = size(SCAPE_data);
    % Create Preview movie/jpg
    if get(handles.makeMIPMovie_checkbox, 'Value') 
        if get(handles.makeMIPMovie_checkbox, 'Value') && (length(ss)==3 || ss(4) == 1)
            filename =  fullfile(handles.directory, [handles.scanName '__HIRES.jpg']);
            if(isdir([handles.directory,'/tiff_stacks']) == 0)
                mkdir([handles.directory,'/tiff_stacks']);
            end
            imgToSave_topMIP = fullfile(handles.directory,'/tiff_stacks/', [handles.scanName '_topMIP_',pf,'_quickMIP.tiff']);
            imgToSave_sideMIP = fullfile(handles.directory,'/tiff_stacks/', [handles.scanName '_sideMIP_', pf,'_quickMIP.tiff']);
                
            tempy = squeeze(max(SCAPE_data(:,3:end-3,:),[],2))';
            imwrite(tempy, imgToSave_topMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');             
            tempy = cat(1,tempy,1000*ones(2,ss(1)));
            tempy = cat(1,tempy, (squeeze(max(SCAPE_data(:,3:end-3,:),[],3))'));
            tmp = max(SCAPE_data(:,3:end-3,:),[],3);

            imwrite(tmp', imgToSave_sideMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');

            [a, b] = hist(reshape(double(tmp),[1,prod(size(tmp))]),200);
            cc1 = b(min(find(a>0.01*mean(a))));
            cc2 = b(max(find(a>0.01*mean(a))));
            tempy = uint8((2^8/(cc2-cc1))*(tempy-cc1));
            imwrite(tempy,filename,'jpg')
            
        else
            if get(handles.incdiff,'Value')==1
                diffy =2; stt = 7;
            else
                diffy = 1;
                stt = 1;
            end
            
            for pp = 1:diffy
                if pp==1
                    filename =  fullfile(handles.directory, [handles.scanName '__dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quikmovie.avi']);
                    imgToSave_topMIP = fullfile(handles.directory,'/tiff_stacks/', [handles.scanName '_topMIP_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP.tiff']);
                    imgToSave_sideMIP = fullfile(handles.directory,'/tiff_stacks/', [handles.scanName '_sideMIP_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP.tiff']);
                
                else
                    filename =  fullfile(handles.directory, [handles.scanName '__DIFF_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'quikmovie.avi']);
                    imgToSave_topMIP = fullfile(handles.directory,'/tiff_stacks/', [handles.scanName '_topMIP_DIFF_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP.tiff']);
                    imgToSave_sideMIP = fullfile(handles.directory,'/tiff_stacks/', [handles.scanName '_sideMIP_DIFF_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP.tiff']);
                
                end
                if exist(filename) 
                    disp('Movie already exists - skipping');
                elseif ~get(handles.makeMIPMovie_checkbox, 'Value') 
                    disp('Movie step skipped!');
                else
                    clear tempy M
                    if pp==1
                        tmp = max(max(SCAPE_data(:,:,:,1),[],4),[],3);
                    else
                        tmp = max(max(max(SCAPE_data(:,:,:,1:end/5),[],4)-min(SCAPE_data(:,:,:,1:end/5),[],4),[],4),[],3);
                    end
                    [a b] = hist(reshape(double(tmp),[1,prod(size(tmp))]),200);
                    cc1 = b(min(find(a>0.01*mean(a))))*.95;
                    cc2 = b(max(find(a>0.01*mean(a))))*.95;
                    if get(handles.makeMIPMovie_checkbox, 'Value')
                        disp('Writing MIP movie..')
                        vidObj = VideoWriter(filename, 'Uncompressed AVI');
                        vidObj.FrameRate = handles.info.info.daq.scanRate/dsf;
                        open(vidObj);
                        map = gray(256);
                    end
                    for volumeCounter = stt:ss(4)
                        numst = 2; % stitch every nth volume
                        clear tempy
                        if pp==1
                            tempy = squeeze(max(SCAPE_data(:,10:end-10,:,volumeCounter),[],2))';
                            imwrite(uint16(tempy), imgToSave_topMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                            
                            tempy = cat(1,tempy,1000*ones(2,ss(1)));
                           
                            tempy = cat(1,tempy, (squeeze(max(SCAPE_data(:,10:end-10,:,volumeCounter),[],3))'));
                            imwrite(uint16(squeeze(max(SCAPE_data(:,3:end-3,:,volumeCounter),[],3))'), imgToSave_sideMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');

                            % disp(volumeCounter)
                            
                        else
                            tempy = squeeze(max(SCAPE_data(:,10:end-10,:,volumeCounter)-SCAPE_data(:,10:end-10,:,volumeCounter-6),[],2))';
                            imwrite(uint16(tempy), imgToSave_topMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                            
                            tempy = cat(1,tempy,1000*ones(2,ss(1)));
                            tempy = cat(1,tempy, (squeeze(max(SCAPE_data(:,10:end-10,:,volumeCounter)-SCAPE_data(:,10:end-10,:,volumeCounter-6),[],3))'));
                            imwrite(uint16(squeeze(max(SCAPE_data(:,3:end-3,:,volumeCounter),[],3))'), imgToSave_sideMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');                        
                        
                        end
                        if get(handles.makeMIPMovie_checkbox, 'Value')
                            M(volumeCounter-stt+1) = im2frame(uint8((256/(cc2-cc1))*(tempy-cc1)),map);
                        end
                    end
                    if get(handles.makeMIPMovie_checkbox, 'Value')
                        writeVideo(vidObj,M);
                        close(vidObj)
                        close(h)
                    end
               
                    disp('done')
                end
            end
        end
    end
    
    if isfield(handles,'tiff')
        if handles.tiff == 1
            if(isdir([handles.directory,'/tiff_stacks']) == 0)
                mkdir([handles.directory,'/tiff_stacks']);
            end
            
            if(isdir([handles.directory,'/tiff_stacks/',handles.scanName]) == 0)
                mkdir([handles.directory,'/tiff_stacks/',handles.scanName]);
            end
            
            filePath = [handles.directory,'/tiff_stacks/',handles.scanName,'/'];
            %filePath = [handles.directory,'/tiff_stacks/'];

            tic
            disp('Start Saving Scan');
            ss = size(SCAPE_data);
            handles.ss = ss;
            conversionFactors = [handles.info.info.cal.ylat; handles.info.info.cal.zdep; handles.info.info.cal.xwid];
            if numel(ss)== 3
                if ~get(handles.orient_lat, 'Value')
                imgToSave = [filePath, 'unCorrected_' handles.scanName '.tiff'];
                fileCounter = 1;
                while (2 == exist(imgToSave, 'file'))
                    imgToSave = [filePath, 'unCorrected_' handles.scanName '_' num2str(fileCounter) '.tiff'];
                    fileCounter = fileCounter + 1;
                end
                else
                imgToSave = [filePath, 'unSkewed_' handles.scanName '.tiff'];
                fileCounter = 1;
                while (2 == exist(imgToSave, 'file'))
                    imgToSave = [filePath, 'unSkewed_' handles.scanName '_' num2str(fileCounter) '.tiff'];
                    fileCounter = fileCounter + 1;
                end
                end
              
                
                if get(handles.checkbox2,'Value')
                    
                    % Coordinate System Correction
                    SCAPE_data = flip(flip(SCAPE_data(:,:,:),2),3);
                    SCAPE_data = permute(SCAPE_data, [3 1 2]);

                    R = imref3d(size(SCAPE_data), conversionFactors(1), conversionFactors(3), conversionFactors(2));
                    
                    %Correct for skew
                    delta = str2num(get(handles.SkewAng,'String'));
                    
                    affineMatrix = [1 0 0 0;
                        0 1 0 0;
                        0 cotd(delta) 1 0;
                        0 0 0 1];
                    tform = affine3d(affineMatrix);
                    [SCAPE_data, ~] = imwarp(SCAPE_data, R, tform);
                    imgToSave = [filePath, handles.scanName '.tiff'];
                    disp('Skew correction applied')
                elseif get(handles.orient_lat, 'Value')
                     % Coordinate System Correction
                    SCAPE_data = flip(flip(SCAPE_data(:,:,:),2),3);
                    SCAPE_data = permute(SCAPE_data, [3 1 2]);
                end

                for j = 1:size(SCAPE_data, 3)
                    imwrite(uint16(SCAPE_data(:, 2:end-1, j))-105, imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                end
            else
                parfor_progress(ss(4));
                
                parfor i  = 1:ss(4)

                    if ~get(handles.orient_lat, 'Value')
                    imgToSave = [filePath, 'unCorrected_' handles.scanName '_t' num2str(i) '.tiff'];
                    fileCounter = 1;
                    while (2 == exist(imgToSave, 'file'))
                        imgToSave = [filePath, 'unCorrected_' handles.scanName '_' num2str(fileCounter) '_t' num2str(i) '.tiff'];
                        fileCounter = fileCounter + 1;
                    end
                    else
                    imgToSave = [filePath, 'unSkewed_' handles.scanName '_t' num2str(i) '.tiff'];
                    fileCounter = 1;
                    while (2 == exist(imgToSave, 'file'))
                        imgToSave = [filePath, 'unSkewed_' handles.scanName '_' num2str(fileCounter) '_t' num2str(i) '.tiff'];
                        fileCounter = fileCounter + 1;
                    end 
                    end
                    
                    temp = SCAPE_data(:,:,:,i);
                    if get(handles.checkbox2,'Value')
                        % Coordinate System Correction
                        temp = flip(flip(temp,2),3);
                        temp = permute(temp, [3 1 2]);
                        R = imref3d(size(temp), conversionFactors(1), conversionFactors(3), conversionFactors(2));
                        
                        % Correct for skew
                        %                     figure(1)
                        %                     imagesc([0:size(SCAPE_data, 1)-1]*conversionFactors(3), ...
                        %                         [0:size(SCAPE_data, 3)-1]*conversionFactors(2), squeeze(max(SCAPE_data(:, 300:350, :), [], 2))')
                        %                     axis image
                        delta = str2num(get(handles.SkewAng,'String'));
                        %delta = pi/2;
                        affineMatrix = [1 0 0 0;
                            0 1 0 0;
                            0 cotd(delta) 1 0;
                            0 0 0 1];
                        tform = affine3d(affineMatrix);
                        [temp, ~] = imwarp(temp, R, tform);
                        % temp = log(double(temp)+1)*2^16/log(2^16)
                        imgToSave = [filePath, handles.scanName '_' num2str(i) '.tiff'];
                    elseif get(handles.orient_lat, 'Value')
                        % Coordinate System Correction
                        temp = flip(flip(temp(:,:,:),2),3);
                        temp = permute(temp, [3 1 2]);
                    end
                    for j = 1:size(temp, 3)
                        imwrite(uint16(temp(:, 1:end, j))-105, imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                    end
                    parfor_progress;
                end
                parfor_progress(0);
                
            end
            disp('Finished Saving Scan');
            toc
            if get(handles.OpeninFiji,'Value')
                addpath('/local_mount/space/enterprise/4/fiji-linux64-20170530/Fiji.app/scripts/')
                javaaddpath '/usr/local/MATLAB/R2016a/java/mij.jar'
                javaaddpath '/usr/local/MATLAB/R2016a/java/ij.jar'
                Miji;
                
                if numel(ss) == 3
                    numtiffs = 1;
                    if get(handles.checkbox2, 'Value')
                        imgToSaveA = [filePath, handles.scanName '.tiff'];
                    else
                        imgToSaveA = [filePath, 'unCorrected_' handles.scanName '.tiff'];
                    end
                else
                    numtiffs = ss(4);
                    if get(handles.checkbox2, 'Value')
                        imgToSaveA = [filePath, handles.scanName '_t1.tiff'];
                    else
                        imgToSaveA = [filePath, 'unCorrected_' handles.scanName '_t1.tiff'];
                    end
                end
                
                info = imfinfo(imgToSaveA);
                z = length(info);
                
                % command1 = strcat('open=[',directory,'] number=', num2str(numtiffs), ' file=.tiff sort');
                % MIJ.run('Image Sequence...', command1)
                
                command1 = strcat('open=[',imgToSaveA,'] number=', num2str(numtiffs), ' starting=1 increment=1 scale=100 file=.tiff');
                MIJ.run('Hypervolume Opener', command1)
                
                command2 = strcat('order=xyczt(default) channels=1 slices=', num2str(z), ' frames=', num2str(numtiffs), ' display=Color');
                MIJ.run('Stack to Hyperstack...', command2);
                
            end
            
        end
    end
    
    [handles] = listbox1_Callback(hObject, eventdata, handles);
    handles.nooutput = 0; 
    guidata(hObject,handles);
    clear SCAPE_data;
    
else
    disp(sprintf('NOTHING TO LOAD: %s',handles.info.info.scanName));
end
