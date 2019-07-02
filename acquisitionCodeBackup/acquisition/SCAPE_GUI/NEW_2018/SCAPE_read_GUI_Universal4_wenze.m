function varargout = SCAPE_read_GUI_Universal4_wenze(varargin)
% SCAPE_READ_GUI_HODOR_CURRENT_PAR MATLAB code for SCAPE_read_GUI_Hodor_Current_par.fig
%      SCAPE_READ_GUI_HODOR_CURRENT_PAR, by itself, creates a new SCAPE_READ_GUI_HODOR_CURRENT_PAR or raises the existing_tiff
%      singleton*.
%
%      H = SCAPE_READ_GUI_HODOR_CURRENT_PAR returns the handle to a new SCAPE_READ_GUI_HODOR_CURRENT_PAR or the handle to
%      the existing_tiff singleton*.
%
%      SCAPE_READ_GUI_HODOR_CURRENT_PAR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCAPE_READ_GUI_HODOR_CURRENT_PAR.M with the given input arguments.
%
%      SCAPE_READ_GUI_HODOR_CURRENT_PAR('Property','Value',...) creates a new SCAPE_READ_GUI_HODOR_CURRENT_PAR or raises the
%      existing_tiff singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SCAPE_read_GUI_Hodor_Current_par_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SCAPE_read_GUI_Hodor_Current_par_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SCAPE_read_GUI_Hodor_Current_par

% Last Modified by GUIDE v2.5 25-May-2017 22:20:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SCAPE_read_GUI_Hodor_Current_par_OpeningFcn, ...
    'gui_OutputFcn',  @SCAPE_read_GUI_Hodor_Current_par_OutputFcn, ...
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


% --- Executes just before SCAPE_read_GUI_Hodor_Current_par is made visible.
function SCAPE_read_GUI_Hodor_Current_par_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SCAPE_read_GUI_Hodor_Current_par (see VARARGIN)
% set(handles.foldername, 'String', 'F:/');
% handles.initialfolder = '/local_mount/space/enterprise/4/SCAPE_data_backup/';
handles.initialfolder = 'D:\';
set(handles.foldername, 'String', handles.initialfolder);
%
% disp('Creating spool file names...')
handles.nooutput = 0;
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
    tempName = tempName(find(tempName ~= ' '));
    tempName = [tempName 'spool.dat'];
    namesOut{spoolCounter} = tempName;
end
namesOut = [{'0000000000spool.dat'} namesOut];
handles.filesToLoad = namesOut;
clear namesOut tempName temp a spoolCounter
guidata(hObject,handles);
% Choose default command line output for SCAPE_read_GUI_Hodor_Current_par
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SCAPE_read_GUI_Hodor_Current_par wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SCAPE_read_GUI_Hodor_Current_par_OutputFcn(hObject, eventdata, handles)
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
    if filesInDataDirectory(i).isdir ==1 & strcmp(filesInDataDirectory(i).name,'movies')==0;
        nameCounter=nameCounter+1;
        runNames{nameCounter} = filesInDataDirectory(i).name;
    end
end

% [SELECTION,OK] = listdlg('ListString',runNames,'ListSize',[400 300]);
handles.directory = directory;
handles.runNames = runNames;
set(handles.listbox1,'Value',1);
set(handles.listbox1,'String',runNames);
[handles] = listbox1_Callback(hObject, eventdata, handles)
guidata(hObject,handles);

% --- Executes on selection change in listbox1.
function [handles] = listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

directory = handles.directory;
filesInDataDirectory = dir(directory);
nameCounter=0;
for i = 1:length(filesInDataDirectory)
    if filesInDataDirectory(i).isdir ==1 && strcmp(filesInDataDirectory(i).name,'movies')==0;
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
%         handles.info.info.daq.stim = rmfield(handles.info.info.daq.stim,'Vout');
    end
    % Creates file path of zyla metadata file (.ini)
    zylaInfoFilePath = fullfile(handles.directory, scanName, 'acquisitionmetadata.ini');
    %     keyboard
    if (2 == exist(zylaInfoFilePath))
        FID = fopen(zylaInfoFilePath, 'r');
        zylaMetaData = fread(FID, '*char')';
        fclose(FID);
        % Reads the number of frames in a spool file
        handles.info.zyla.numDepths = str2double(zylaMetaData(12+(strfind(zylaMetaData, 'AOIHeight = ')):(strfind(zylaMetaData, 'AOIWidth = ')-1)));
        handles.info.zyla.numLatPix = str2double(zylaMetaData((11+strfind(zylaMetaData, 'AOIWidth = ')):(strfind(zylaMetaData, 'AOIStride = ')-1)));
        handles.info.zyla.imageBytes = str2double(zylaMetaData((17+strfind(zylaMetaData, 'ImageSizeBytes = ')):(strfind(zylaMetaData, '[multiimage]')-1)));
        handles.info.zyla.numFramesPerSpool = str2double(zylaMetaData((16+strfind(zylaMetaData, 'ImagesPerFile = ')):end));
        startIndex = strfind(zylaMetaData, 'ImageSizeBytes') + length('ImageSizeBytes = ');
        handles.info.zyla.ImageSize = str2double(zylaMetaData(startIndex : startIndex + 7));
    else
        disp('No metadata file found for this run');
    end

    %generate string of info:
    clear infotext
    infotext{1} = sprintf('Scan Rate: %.2f VPS',handles.info.info.daq.scanRate);
    infotext{2} = sprintf('Scan Duration: %.2f secs',handles.info.info.daq.scanLength);
    infotext{3} = sprintf('X-steps: %d, X-range: %.4f V',handles.info.info.daq.pixelsPerLine,handles.info.info.daq.scanAngle);
    infotext{4} = sprintf('Camera Frame Rate: %.2f Hz',handles.info.info.daq.pixelrate);
    infotext{5} = sprintf('Camera height: %d, Width: %d',handles.info.zyla.numDepths,handles.info.zyla.numLatPix);
    
    if(isfield(handles.info.info,'sawtooth'));
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
    scandate = datenum(handles.info.info.scanStartTimeApprox);
    calpath = 'C:\Users\Goat\Documents\MATLAB\SCAPE\Load Code\';
    if exist([calpath 'SCAPE_calibrations.txt'],'file')==2
        fid = fopen([calpath 'SCAPE_calibrations.txt'],'r');
        fgetl(fid);
        %         fgetl(fid);
        m=1;
        a{m} = fgetl(fid);
        while a{m}~=-1
            m=m+1;
            a{m}= fgetl(fid);
            p = strfind(a{m},'	');
            try sdate = datenum(a{m}(1:p(1)-1));
                edate = datenum(a{m}(p(1)+1:p(2)-1));
                if scandate>=sdate && scandate<=edate
                    handles.info.info.cal = [];
                    handles.info.info.cal.ylat = str2num(a{m}(p(2)+1:p(3)-1));
                    handles.info.info.cal.zdep = str2num(a{m}(p(3)+1:p(4)-1));
                    handles.info.info.cal.xfact = str2num(a{m}(p(4)+1:end));
                    handles.info.info.cal.xwid = handles.info.info.cal.xfact*handles.info.info.daq.scanAngle/handles.info.info.daq.pixelsPerLine;
                    a{m} = -1;
                    success = 1;
                end
            catch
                success = 0;
            end
        end
        fclose(fid)
        if success ==1
            infotext{8} = sprintf('Calibration factors: y: %.2f  x: %.1f z: %.1f um/pix',handles.info.info.cal.ylat, handles.info.info.cal.xwid,handles.info.info.cal.zdep);
            infotext{9} = sprintf('Field of view: y: %.2f  x: %.1f z: %.1f um',handles.info.info.cal.ylat*handles.info.zyla.numLatPix, handles.info.info.cal.xwid*handles.info.info.daq.pixelsPerLine,handles.info.zyla.numDepths*handles.info.info.cal.zdep);
        else
            disp('! calibration factors for this DATE not found');
            handles.info.info.cal.ylat = 1.39;
            handles.info.info.cal.zdep = 	1.14;
            handles.info.info.cal.xfact = 	323.408;
            handles.info.info.cal.xwid = handles.info.info.cal.xfact*handles.info.info.daq.scanAngle/handles.info.info.daq.pixelsPerLine;
            handles.info.info.cal2 = sprintf('Get recent vals! Using xwid=%0.2f ylat=1.24 zdep=1.22', handles.info.info.cal.xwid);
            infotext{8} = sprintf('Calibration: %s',handles.info.info.cal2);
        end
    else
        disp('! calibration factor FILE not found');
        handles.info.info.cal.ylat = 1.39;
        handles.info.info.cal.zdep = 	1.14;
        handles.info.info.cal.xfact = 	323.408;
        handles.info.info.cal.xwid = handles.info.info.cal.xfact*handles.info.info.daq.scanAngle/handles.info.info.daq.pixelsPerLine;
        handles.info.info.cal2 = sprintf('Get recent vals! Using xwid=%0.2f ylat=1.24 zdep=1.22', handles.info.info.cal.xwid);
        infotext{8} = sprintf('Calibration: %s',handles.info.info.cal2);
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
    
    
    
    m=0;
    existing{1} = 'no TIFF stacks found';
    m=m+1;
    %  infotext{9} = 'existing_tiff outputs found:';
    if isdir([handles.directory,'/tiff_stacks'])
        tiffs = dir([handles.directory,'/tiff_stacks']);
        for i = 1:length(tiffs)
            if strfind(tiffs(i).name,scanName);
                numtiffs = dir([handles.directory,'/tiff_stacks/',tiffs(i).name]);
                existing{m} = sprintf('%s  (%d files)',tiffs(i).name, length(numtiffs)-2);
                m=m+1;
            end
        end
    end
    set(handles.existing_tiff,'String',existing);
    
    m=1;
    set(handles.existing_avi,'Value',1);
    existing_avi{1} = 'none';
    avis = dir([handles.directory,'/*.avi']);
    for i = 1:length(avis)
        if strfind(avis(i).name,scanName);
            existing_avi{m} = sprintf('%s',avis(i).name);
            m=m+1;
        end
    end
    set(handles.existing_avi,'String',existing_avi);
    %
    %  if exist(filename)
    %       infotext(9} = sprintf('Tiff files found: %s', tifnames);
    %  end
    set(handles.text2,'String',infotext);
    
    directory = handles.directory;
    filesInDataDirectory = dir(directory);
    nameCounter=0;
    for i = 1:length(filesInDataDirectory)
        if filesInDataDirectory(i).isdir ==1 & strcmp(filesInDataDirectory(i).name,'movies')==0;
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
function existing_avi_Callback(hObject, eventdata, handles)
% hObject    handle to existing_avi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns existing_avi contents as cell array
%        contents{get(hObject,'Value')} returns selected item from existing_avi
vids = get(handles.existing_avi,'String');
vidchose = get(handles.existing_avi,'Value');
if strcmp(char(vids(vidchose)),'none')==0
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
    % for i = 1:100
    %     imagesc()
    % end
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

% dataDirectory = get(handles.foldername,'String');
% if(isdir([dataDirectory,'/tiff_stacks']) == 0)
%     mkdir([dataDirectory,'/tiff_stacks']);
% end
%
% if(isdir([dataDirectory,'/tiff_stacks/',scanName]) == 0)
%     mkdir([dataDirectory,'/tiff_stacks/',scanName]);
% end
%
% filePath = [dataDirectory,'/tiff_stacks/',scanName,'/'];
% tic
% h = waitbar(0, 'Please Wait.');
% for i = 1:handles.ss(4);%scans
%         m=m+1;
%         imgToSave = [filepath,'RGB_' scanName, '_', num2str(i) '.tiff'];
%         %normalize each volume to its max
%           [red green] = go_Callback(hObject, eventdata, handles,i);
%
%         for j = 1:handles.ss(3)
%             clear rgb
%             rgb(:,:,1) = squeeze(scr*(double(red(:,:,j,i))-minr));
%             rgb(:,:,2) = squeeze(scg*(double(green(:,:,j,i))-ming));
%             rgb(:,:,3) = zeros(size(rgb(:,:,1)));
%             % temp1 = squeeze(handles.inputdata(:,:,:,i));
%             temp11 = squeeze(uint8(256*rgb));%(2^16)*temp1/(max(max(max(temp1))));
%             imwrite(temp11(1:crp,:,:), imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
%         end
% for i = 1:1:size(moviesub, 3)
%     %imwrite(uint16(1000*log10(double(squeeze(skewedMoviesub_il(:, i, :))))), imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
%     imwrite(squeeze(moviesub(200:end-200, :, i)), imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
%     %imwrite(squeeze(abs(X(:, i, :))), imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
%     waitbar(i/size(moviesub, 2));
% end
% close (h)
% disp('Finished Saving Scan');
% toc

% % Wenze Ver
% loaddata_Callback(hObject, eventdata, handles)
% [handles] = listbox1_Callback(hObject, eventdata, handles)
% dataDirectory = get(handles.foldername,'String');
% if(isdir([dataDirectory,'/tiff_stacks']) == 0)
%     mkdir([dataDirectory,'/tiff_stacks']);
% end
%
% if(isdir([dataDirectory,'/tiff_stacks/',handles.scanName]) == 0)
%     mkdir([dataDirectory,'/tiff_stacks/',handles.scanName]);
% end
%
% filePath = [dataDirectory,'/tiff_stacks/',handles.scanName,'/'];
%
% tic
% h = waitbar(0, 'Please Wait.');
% ss = size(SCAPE_data);
% if numel(ss)==3; ss(4) = 1; end
% for i = 1:ss(4) %scans
%     imgToSave = [filePath, 'unCorrected_' handles.scanName '_' num2str(i) '.tiff'];
%     for j = 1:size(SCAPE_data, 3)
%         imwrite(uint16(SCAPE_data(:, 1:end, j, i)), imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
%     end
%     waitbar(i/ss(4));
%
% end
% close (h)
% disp('Finished Saving Scan');
% toc
% dataDirectory = get(handles.foldername,'String');
% if(isdir([dataDirectory,'/tiff_stacks']) == 0)
%     mkdir([dataDirectory,'/tiff_stacks']);
% end
%
% if(isdir([dataDirectory,'/tiff_stacks/',scanName]) == 0)
%     mkdir([dataDirectory,'/tiff_stacks/',scanName]);
% end


handles.nooutput = 1;
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
            [handles] = listbox1_Callback(hObject, eventdata, handles)
            handles.tiff = 1;
            guidata(hObject,handles);
            filename =  fullfile(handles.directory, [handles.scanName '__dsf',num2str(get(handles.dsf,'String')), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',handles.pf,'_quikmovie.avi']);
            %             if exist(filename)
            %                 disp(sprintf('Movie already exists - skipping (%s)',handles.scanName));
            %             else
            loaddata_Callback(hObject, eventdata, handles)
            [handles] = listbox1_Callback(hObject, eventdata, handles)
        end
    end
    handles.nooutput = 0;
    handles.tiff = 0;
    guidata(hObject,handles);
end

% --- Executes on button press in makeavi.
function makeavi_Callback(hObject, eventdata, handles)
% hObject    handle to makeavi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nooutput = 1;
if strcmp(get(handles.secs,'String'),'all')==0
    %     [button] = questdlg('Total duration is not set to ALL - continue?','Yes','Woops - ok, stop');
    disp('WARNING: Total duration is not set to ALL ');
end
button = 'Yes';
handles.tiff = 0;
if strcmp(button,'Yes')
    [SELECTION,OK] = listdlg('ListString',handles.runNames,'ListSize',[400 300],'InitialValue',get(handles.listbox1,'Value'));
    for i = 1:length(SELECTION)
        if strcmp(handles.runNames(SELECTION(i)),'.')==0 && strcmp(handles.runNames(SELECTION(i)),'..')==0
            set(handles.listbox1,'Value',SELECTION(i));
            [handles] = listbox1_Callback(hObject, eventdata, handles)
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
[handles] = listbox1_Callback(hObject, eventdata, handles)
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

numFramesPerSpool = handles.info.zyla.numFramesPerSpool;
numDepths = handles.info.zyla.numDepths;
numLatPix = handles.info.zyla.numLatPix;


% numSpoolfilesPerVolume = ceil(handles.info.info.daq.pixelsPerLine/handles.info.zyla.numFramesPerSpool);
numSpoolfilesPerVolume = ceil(framesPerScan/numFramesPerSpool);
numSpoolfilesPerVolume1 = (framesPerScan/numFramesPerSpool);

numActual_spoolfiles = length(dir([handles.directory,'/',handles.scanName,'/*.dat']));
if numSpoolfilesPerVolume1*scanRate*numSecondsToLoad>numActual_spoolfiles
    numSecondsToLoad = floor(numActual_spoolfiles/(numSpoolfilesPerVolume1*scanRate));
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
    title('CROPPING ONLY Y DIMENSION')
    clear framey;
    [YY xx] = ginput(2);
    Ycrop(1) = ceil(min(YY));
    Ycrop(2) = floor(max(YY));
    %     cropping = 1
else
    %     cropping = 0;
    Ycrop = [1 numRows-10];
    xx = [1 handles.info.info.daq.pixelsPerLine];
end
handles.info.info.Ycrop = Ycrop;

% pre-loaded spool files (add check that there are enough)
filesToLoad = handles.filesToLoad;

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
    startFrameInVolumeToLoad = frame_starts(floor(volumesToLoad));
    
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
    
    if get(handles.makeDouble_checkbox, 'Value');
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
%     parfor_progress(numFileToLoad);
%     keyboard
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
            SCAPE_data(:,:,:,spoolFileCounter) = double(rawData(Ycrop(1):Ycrop(2),1:end-3,:));
        else
            SCAPE_data(:,:,:,spoolFileCounter) = rawData(Ycrop(1):Ycrop(2),1:end-3,:);
        end
%         parfor_progress;
        disp(spoolFileCounter);
    end
%     parfor_progress(0);
% keyboard
    % % Venkata: The need for this step cuts the amount of data we can load
    % into memory in half. Ridiculous...
    
    if volumesToLoad == 1
        SCAPE_data = reshape(SCAPE_data,[[Ycrop(2)-Ycrop(1)]+1, numColumns-3, numFramesPerSpool*numFileToLoad]);
        SCAPE_data = reshape(SCAPE_data(:,:,index),[[Ycrop(2)-Ycrop(1)]+1, numColumns-3, framesPerScan-2, numVolumeToLoad]);
    else
        SCAPE_data = reshape(SCAPE_data,[[Ycrop(2)-Ycrop(1)]+1, numColumns-3, numFramesPerSpool*numFileToLoad]);
        SCAPE_data = reshape(SCAPE_data(:,:,index),[[Ycrop(2)-Ycrop(1)]+1, numColumns-3, framesPerScan-2, numVolumeToLoad]);
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
        assignin('base','info',handles.info.info);
        assignin('base','zyla',handles.info.zyla);
        h = waitbar(1,'SCAPE data written to workspace');
        pause(2);
        close(h)
    else
        h = waitbar(1,'Writing movie');
        pause(2);
        close(h)
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
        if length(ss)==3 || ss(4) == 1
            filename =  fullfile(handles.directory, [handles.scanName '__HIRES.jpg']);
            tempy = squeeze(max(SCAPE_data(:,3:end-3,:),[],2))';
            tempy = cat(1,tempy,1000*ones(2,ss(1)));
            tempy = cat(1,tempy, (squeeze(max(SCAPE_data(:,3:end-3,:),[],3))'));
            tmp = max(SCAPE_data(:,3:end-3,:),[],3);
            [a, b] = hist(reshape(double(tmp),[1,prod(size(tmp))]),200);
            cc1 = b(min(find(a>0.01*mean(a))));
            cc2 = b(max(find(a>0.01*mean(a))));
            tempy = uint8((2^8/(cc2-cc1))*(tempy-cc1));
            imwrite(tempy,filename,'jpg')
        else
            
            if get(handles.incdiff,'Value')==1;
                diffy =2; stt = 7;
            else
                diffy = 1;
                stt = 1;
            end
            
            for pp = 1:diffy;
                dataDirectory = get(handles.foldername,'String');
                if(isdir([dataDirectory,'/tiff_stacks']) == 0)
                    [dataDirectory,'/tiff_stacks']
                    mkdir([dataDirectory,'/tiff_stacks']);
                end
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
                    %                 close(h);
                elseif (get(handles.makeMIPMovie_checkbox, 'Value') == 0)
                    disp('Movie step skipped!');
                else
                    disp('Writing MIP movie..')
                    clear tempy M
                    if pp==1
                        tmp = max(max(SCAPE_data(:,:,:,1),[],4),[],3);
                    else
                        tmp = max(max(max(SCAPE_data(:,:,:,1:end/5),[],4)-min(SCAPE_data(:,:,:,1:end/5),[],4),[],4),[],3);
                    end
                    [a b] = hist(reshape(double(tmp),[1,prod(size(tmp))]),200);
                    cc1 = b(min(find(a>0.01*mean(a))))*.95;
                    cc2 = b(max(find(a>0.01*mean(a))))*.95;
                    vidObj = VideoWriter(filename, 'Uncompressed AVI');
                    vidObj.FrameRate = handles.info.info.daq.scanRate/dsf;
                    open(vidObj);
                    map = gray(256);
                    
                    for volumeCounter = stt:ss(4)
                        clear tempy
                        if pp==1
                            tempy = squeeze(max(SCAPE_data(:,3:end-3,:,volumeCounter),[],2))';
                            imwrite(uint16(tempy), imgToSave_topMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                            tempy = cat(1,tempy,1000*ones(2,ss(1)));                            
                            tempy = cat(1,tempy, (squeeze(max(SCAPE_data(:,3:end-3,:,volumeCounter),[],3))'));
                            imwrite(uint16(squeeze(max(SCAPE_data(:,3:end-3,:,volumeCounter),[],3))'), imgToSave_sideMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');

                            %             disp(volumeCounter)
                        else
                            tempy = squeeze(max(SCAPE_data(:,3:end-3,:,volumeCounter)-SCAPE_data(:,3:end-3,:,volumeCounter-6),[],2))';
                            imwrite(uint16(tempy), imgToSave_topMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                            tempy = cat(1,tempy,1000*ones(2,ss(1)));
                            tempy = cat(1,tempy, (squeeze(max(SCAPE_data(:,3:end-3,:,volumeCounter)-SCAPE_data(:,3:end-3,:,volumeCounter-6),[],3))'));
                            imwrite(uint16(squeeze(max(SCAPE_data(:,3:end-3,:,volumeCounter),[],3))'), imgToSave_sideMIP, 'tif', 'Compression', 'none', 'WriteMode', 'Append');                        
                        end
                        
                        M(volumeCounter-stt+1) = im2frame(uint8((256/(cc2-cc1))*(tempy-cc1)),map);
                    end
                    writeVideo(vidObj,M);
                    close(vidObj)
                    
                    disp('done')
                end
            end
        end
    end
    
    
    if isfield(handles,'tiff')
        if handles.tiff == 1
            dataDirectory = get(handles.foldername,'String');
            if(isdir([dataDirectory,'/tiff_stacks']) == 0)
                [dataDirectory,'/tiff_stacks']
                mkdir([dataDirectory,'/tiff_stacks']);
            end
            
            if(isdir([dataDirectory,'/tiff_stacks/',handles.scanName]) == 0)
                mkdir([dataDirectory,'/tiff_stacks/',handles.scanName]);
            end
            
            filePath = [dataDirectory,'/tiff_stacks/',handles.scanName,'/'];
            
            tic
            disp('Start Saving Scan');
            ss = size(SCAPE_data);
            conversionFactors = [handles.info.info.cal.ylat; handles.info.info.cal.zdep; handles.info.info.cal.xwid];
            if numel(ss)== 3
                imgToSave = [filePath, 'unCorrected_' handles.scanName '.tiff'];
                fileCounter = 1;
                while (2 == exist(imgToSave, 'file'))
                    imgToSave = [filePath, 'unCorrected_' handles.scanName '_' num2str(fileCounter) '.tiff'];
                    fileCounter = fileCounter + 1;
                end
                
                if get(handles.checkbox2,'Value')
                    
                    % Coordinate System Correction
                    SCAPE_data = flip(flip(SCAPE_data(:,:,:),2),3);
                    SCAPE_data = permute(SCAPE_data, [3 1 2]);
                    R = imref3d(size(SCAPE_data), conversionFactors(3), conversionFactors(1), conversionFactors(2));
                    
                    % Correct for skew
                    
                    %                        figure(99); subplot(1,2,1)
                    %                        imagesc([0:size(SCAPE_data, 1)-1]*conversionFactors(3),[0:size(SCAPE_data, 3)-1]*conversionFactors(2), log(double(squeeze(max(SCAPE_data(:,round(end/3):round(2*end/3), :), [], 2))))')
                    %                        colormap gray; axis image
                    
                    delta = str2num(get(handles.SkewAng,'String'))*pi/180;
                    
                    affineMatrix = [1 0 0 0;
                        0 1 0 0;
                        0 cot(delta) 1 0;
                        0 0 0 1];
                    tform = affine3d(affineMatrix);
                    [SCAPE_data, ~] = imwarp(SCAPE_data, R, tform);
                    
                    %                     figure(99); subplot(1,2,2)
                    %                     imagesc([0:size(SCAPE_data, 1)-1]*conversionFactors(3),[0:size(SCAPE_data, 3)-1]*conversionFactors(2), log(double(squeeze(max(SCAPE_data(:,round(end/3):round(2*end/3), :), [], 2))))')
                    %                     colormap gray; axis image
                    imgToSave = [filePath, handles.scanName '.tiff'];
                    disp('Skew correction applied')
                end
                
                for j = 1:size(SCAPE_data, 3)
                    imwrite(uint16(SCAPE_data(:, 1:end, j)), imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                end
                
            else
            parfor_progress(ss(4));
                
                parfor (i  =  1:ss(4),8) %scans
                    imgToSave = [filePath, 'unCorrected_' handles.scanName '_t' num2str(i) '.tiff'];
                    fileCounter = 1;
                    while (2 == exist(imgToSave, 'file'))
                        imgToSave = [filePath, 'unCorrected_' handles.scanName '_' num2str(fileCounter) '_t' num2str(i) '.tiff'];
                        fileCounter = fileCounter + 1;
                    end
                    temp = SCAPE_data(:,:,:,i);
                    if get(handles.checkbox2,'Value')
                        % Coordinate System Correction
                        temp = flip(flip(temp,2),3);
                        temp = permute(temp, [3 1 2]);
                        R = imref3d(size(temp), conversionFactors(3), conversionFactors(1), conversionFactors(2));
                        
                        % Correct for skew
                        %                     figure(1)
                        %                     imagesc([0:size(SCAPE_data, 1)-1]*conversionFactors(3), ...
                        %                         [0:size(SCAPE_data, 3)-1]*conversionFactors(2), squeeze(max(SCAPE_data(:, 300:350, :), [], 2))')
                        %                     axis image
                        delta = str2num(get(handles.SkewAng,'String'))*pi/180;
                        %delta = pi/2;
                        affineMatrix = [1 0 0 0;
                            0 1 0 0;
                            0 cot(delta) 1 0;
                            0 0 0 1];
                        tform = affine3d(affineMatrix);
                        [temp, ~] = imwarp(temp, R, tform);
                        imgToSave = [filePath, handles.scanName '_' num2str(i) '.tiff'];
                    end
                    for j = 1:size(temp, 3)
                        imwrite(uint16(temp(:, 1:end, j)), imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                    end
               parfor_progress;
                end
              parfor_progress(0);
                
            end
            disp('Finished Saving Scan');
        end
    end
    
    [handles] = listbox1_Callback(hObject, eventdata, handles);
    
    clear SCAPE_data;
    
else
    disp(sprintf('NOTHING TO LOAD: %s',handles.info.info.scanName));
    keyboard
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


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(~, ~, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Call Splitcols
SCAPE_splitcolsV5(evalin('base', 'SCAPE_data'),evalin('base', 'info'))


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
handles.info.info.cal.ylat = 1.24;
handles.info.info.cal.zdep = 	1.22;
handles.info.info.cal.xfact = 	327.28;
conversionFactors = [handles.info.info.cal.ylat; handles.info.info.cal.zdep; handles.info.info.cal.xwid];

SCAPE_data = flip(flip(SCAPE_data(:,:,:),2),3);
SCAPE_data = permute(SCAPE_data, [3 1 2]);
R = imref3d(size(SCAPE_data), conversionFactors(3), conversionFactors(1), conversionFactors(2));

figure(99); subplot(2,1,1)
imagesc([0:size(SCAPE_data, 1)-1]*conversionFactors(3),[0:size(SCAPE_data, 3)-1]*conversionFactors(2), log(double(squeeze(max(SCAPE_data(:,round(end/3):round(2*end/3), :), [], 2))))')
colormap gray; axis image; title('Raw')

delta = str2num(get(handles.SkewAng,'String'))*pi/180;

affineMatrix = [1 0 0 0;
    0 1 0 0;
    0 cot(delta) 1 0;
    0 0 0 1];
tform = affine3d(affineMatrix);
[SCAPE_data, ~] = imwarp(SCAPE_data, R, tform);

figure(99); subplot(2,1,2)
imagesc([0:size(SCAPE_data, 1)-1]*conversionFactors(3),[0:size(SCAPE_data, 3)-1]*conversionFactors(2), log(double(squeeze(max(SCAPE_data(:,round(end/3):round(2*end/3), :), [], 2))))')
colormap gray; axis image; title('Skew Corrected')
imgToSave = [filePath, handles.scanName '.tiff'];
disp('Skew correction applied')


% --- Executes on button press in dual_color_tiff_checkbox.
function dual_color_tiff_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to dual_color_tiff_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dual_color_tiff_checkbox


% --- Executes on button press in dualcolorpreview.
function dualcolorpreview_Callback(hObject, eventdata, handles)
% hObject    handle to dualcolorpreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nooutput = 1;
if strcmp(get(handles.secs,'String'),'all')==0
    %     [button] = questdlg('Total duration is not set to ALL - continue?','Yes','Woops - ok, stop');
    disp('WARNING: Total duration is not set to ALL ');
end
button = 'Yes';
handles.tiff = get(handles.dual_color_tiff_checkbox, 'Value');
if strcmp(button,'Yes')
    [SELECTION,OK] = listdlg('ListString',handles.runNames,'ListSize',[400 300],'InitialValue',get(handles.listbox1,'Value'));
    for i = 1:length(SELECTION)
        if strcmp(handles.runNames(SELECTION(i)),'.')==0 && strcmp(handles.runNames(SELECTION(i)),'..')==0
            set(handles.listbox1,'Value',SELECTION(i));
            [handles] = listbox1_Callback(hObject, eventdata, handles)
            guidata(hObject,handles);
            filename =  fullfile(handles.directory, [handles.scanName '__dsf',num2str(get(handles.dsf,'String')), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',handles.pf,'_dualmovie.avi']);
            if exist(filename)
                disp(sprintf('Movie already exists - skipping (%s)',handles.scanName));
            else
                dual_loaddata_Callback(hObject, eventdata, handles)
                [handles] = listbox1_Callback(hObject, eventdata, handles)
            end
        end
    end
    handles.nooutput = 0;
    guidata(hObject,handles);
end

% dualcolor version load code
function dual_loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to loaddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hh2 = waitbar(1,'Working');
[handles] = listbox1_Callback(hObject, eventdata, handles)
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

numFramesPerSpool = handles.info.zyla.numFramesPerSpool;
numDepths = handles.info.zyla.numDepths;
numLatPix = handles.info.zyla.numLatPix;


% numSpoolfilesPerVolume = ceil(handles.info.info.daq.pixelsPerLine/handles.info.zyla.numFramesPerSpool);
numSpoolfilesPerVolume = ceil(framesPerScan/numFramesPerSpool);
numSpoolfilesPerVolume1 = (framesPerScan/numFramesPerSpool);

numActual_spoolfiles = length(dir([handles.directory,'/',handles.scanName,'/*.dat']));
if numSpoolfilesPerVolume1*scanRate*numSecondsToLoad>numActual_spoolfiles
    numSecondsToLoad = floor(numActual_spoolfiles/(numSpoolfilesPerVolume1*scanRate));
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
    title('CROPPING ONLY Y DIMENSION')
    clear framey;
    [YY xx] = ginput(2);
    Ycrop(1) = ceil(min(YY));
    Ycrop(2) = floor(max(YY));
    %     cropping = 1
else
    %     cropping = 0;
    Ycrop = [1 numRows-10];
    xx = [1 handles.info.info.daq.pixelsPerLine];
end
handles.info.info.Ycrop = Ycrop;

% pre-loaded spool files (add check that there are enough)
filesToLoad = handles.filesToLoad;

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
    
    if get(handles.makeDouble_checkbox, 'Value');
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
            SCAPE_data(:,:,:,spoolFileCounter) = double(rawData(Ycrop(1):Ycrop(2),1:end-3,:));
        else
            SCAPE_data(:,:,:,spoolFileCounter) = rawData(Ycrop(1):Ycrop(2),1:end-3,:);
        end
        parfor_progress;
    end
    parfor_progress(0);
    
    if volumesToLoad == 1
        SCAPE_data = reshape(SCAPE_data,[[Ycrop(2)-Ycrop(1)]+1, numColumns-3, numFramesPerSpool*numFileToLoad]);
        SCAPE_data = reshape(SCAPE_data(:,:,index),[[Ycrop(2)-Ycrop(1)]+1, numColumns-3, framesPerScan-2, numVolumeToLoad]);
    else
        SCAPE_data = reshape(SCAPE_data,[[Ycrop(2)-Ycrop(1)]+1, numColumns-3, numFramesPerSpool*numFileToLoad]);
        SCAPE_data = reshape(SCAPE_data(:,:,index),[[Ycrop(2)-Ycrop(1)]+1, numColumns-3, framesPerScan-2, numVolumeToLoad]);
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
        assignin('base','info',handles.info.info);
        assignin('base','zyla',handles.info.zyla);
        h = waitbar(1,'SCAPE data written to workspace');
        pause(2);
        close(h)
    else
        h = waitbar(1,'Writing movie');
        pause(2);
        close(h)
    end
    
    % Flips even frames to account for bidirectional scan. Also waits until second volume is loaded and flipped to introduce the pixel shift.
    if (saw == 0 && (round(dsf/2)== dsf/2)==0)
        SCAPE_data(:, :, :, 1:2:end) = flipdim(SCAPE_data(:, :, :, 1:2:end), 3);
    end
    
    % close(h)
    clear M
    ss = size(SCAPE_data);
    sn = ss;
    % skew correction
    if get(handles.checkbox2,'Value')

        green_y_project=max(SCAPE_data(1:(end/2),:,:,floor(sn(4)/10)),[],3);
        red_y_project=max(SCAPE_data((end/2)+1:end,:,:,floor(sn(4)/10)),[],3);
        tform_y=imregcorr(green_y_project,red_y_project,'rigid');
        tform_y.T
        angle = str2num(get(handles.SkewAng,'String'));
        absangle = 1/tan(2*pi*angle/360);
        conversionFactors = [handles.info.info.cal.ylat; handles.info.info.cal.zdep; handles.info.info.cal.xwid];

        shps = absangle*conversionFactors(1)/conversionFactors(3);

        skewshift = zeros(sn(1),sn(2),ceil(shps*sn(2))+sn(3),sn(4),'uint16');
%         m=0; 
        disp('Skew Correction Start')
tic
        for i = ([1:sn(2)])
%             m=m+1;
            bba = floor(i*shps);
            temp = double(squeeze(SCAPE_data(:,i,:,:)));
            skewshift(:,i,[1:sn(3)]+bba,:) = uint16(interp3([1:sn(3)],[1:sn(1)]',[1:sn(4)],temp,[1:sn(3)]-i*shps+bba,[1:sn(1)]',[1:sn(4)],'linear'));
%             disp(i)
        end
        disp('Skew Correction Done')
        clear SCAPE_data   
        toc
    else
        skewshift = SCAPE_data;
        clear SCAPE_data 
    end
    % Green and red split
    ss = size(skewshift);
%     zz = tform_y.T(3,1);
%     yy = tform_y.T(3,2);    
    zz = 2.2;
    yy = 19;
    green = zeros(ss(1)/2,ss(2),ss(3),ss(4),'uint16');
    % shps = -0.635
    disp('Dual Color Split Start')
tic
    for i = ([1:ss(3)])
        temp = double(squeeze(skewshift(1:end/2,:,i,:)));
        green(:,:,i,:) = uint16(interp3([1:ss(2)],[1:ss(1)/2]',[1:ss(4)],temp,[1:ss(2)]-zz,-yy+[1:ss(1)/2]',[1:ss(4)],'linear'));
%         disp(i)
    end
    red = skewshift((end/2)+1:end,:,:,:);    
    disp('Dual Color Split Done')
    ss = size(green);
    clear skewshift
    toc
    
    vidlength = ss(4);
    % Create Preview movie/jpg
    if get(handles.makeMIPMovie_checkbox, 'Value')
        if length(ss)==3 || ss(4) == 1
            % for HR
            filename =  fullfile(handles.directory, [handles.scanName '__HIRES.jpg']);
            tempy = squeeze(max(SCAPE_data(:,3:end-3,:),[],2))';
            tempy = cat(1,tempy,1000*ones(2,ss(1)));
            tempy = cat(1,tempy, (squeeze(max(SCAPE_data(:,3:end-3,:),[],3))'));
            tmp = max(SCAPE_data(:,3:end-3,:),[],3);
            [a, b] = hist(reshape(double(tmp),[1,prod(size(tmp))]),200);
            cc1 = b(min(find(a>0.01*mean(a))));
            cc2 = b(max(find(a>0.01*mean(a))));
            tempy = uint8((2^8/(cc2-cc1))*(tempy-cc1));
            imwrite(tempy,filename,'jpg')
        else
            
            if get(handles.incdiff,'Value')==1;
                diffy =2; stt = 7;
            else
                diffy = 1;
                stt = 1;
            end
            
            for pp = 1:diffy;
                dataDirectory = get(handles.foldername,'String');
                if(isdir([dataDirectory,'/dual_color_preview_new']) == 0)
                    [dataDirectory,'/dual_color_preview_new']
                    mkdir([dataDirectory,'/dual_color_preview_new']);
                end
                if pp==1
                    filename =  fullfile(handles.directory,'/dual_color_preview_new', [handles.scanName '__dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quikmovie_dual.avi']);
                    imgToSave_topMIP_G = fullfile(handles.directory,'/dual_color_preview_new/', [handles.scanName '_topMIP_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP_G.tiff']);
                    imgToSave_sideMIP_G = fullfile(handles.directory,'/dual_color_preview_new/', [handles.scanName '_sideMIP_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP_G.tiff']);
                    imgToSave_topMIP_R = fullfile(handles.directory,'/dual_color_preview_new/', [handles.scanName '_topMIP_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP_R.tiff']);
                    imgToSave_sideMIP_R = fullfile(handles.directory,'/dual_color_preview_new/', [handles.scanName '_sideMIP_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP_R.tiff']);
                else
                    filename =  fullfile(handles.directory, '/dual_color_preview_new',[handles.scanName '__log_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'quikmovie_dual.avi']);
                    imgToSave_topMIP_G = fullfile(handles.directory,'/dual_color_preview_new/', [handles.scanName '_topMIP_log_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP_G.tiff']);
                    imgToSave_sideMIP_G = fullfile(handles.directory,'/dual_color_preview_new/', [handles.scanName '_sideMIP_log_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP_G.tiff']);
                    imgToSave_topMIP_R = fullfile(handles.directory,'/dual_color_preview_new/', [handles.scanName '_topMIP_log_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP_R.tiff']);
                    imgToSave_sideMIP_R = fullfile(handles.directory,'/dual_color_preview_new/', [handles.scanName '_sideMIP_log_dsf',num2str(dsf), '_',num2str(handles.info.info.daq.scanRate),'VPS_',get(handles.secs,'String'),'secs_',pf,'_quickMIP_R.tiff']);
                end
                
                if exist(filename)
                    disp('Movie already exists - skipping');
                    %                 close(h);
                elseif (get(handles.makeMIPMovie_checkbox, 'Value') == 0)
                    disp('Movie step skipped!');
                else
                    disp('Writing MIP movie..')
                    clear tempy M
%                     if pp==1
%                         tmp = max(max(SCAPE_data(:,:,:,1),[],4),[],3);
%                     else
%                         tmp = max(max(max(SCAPE_data(:,:,:,1:end/5),[],4)-min(SCAPE_data(:,:,:,1:end/5),[],4),[],4),[],3);
%                     end
%                     [a b] = hist(reshape(double(tmp),[1,prod(size(tmp))]),200);
%                     cc1 = b(min(find(a>0.01*mean(a))))*.95;
%                     cc2 = b(max(find(a>0.01*mean(a))))*.95;
                    vidObj = VideoWriter(filename, 'Uncompressed AVI');
                    vidObj.FrameRate = handles.info.info.daq.scanRate/dsf;
                    open(vidObj);
%                     map = gray(256);
                    
                    for volumeCounter = stt:vidlength
                        clear tempy rgb rgb2
                        if pp==1
%                             disp(volumeCounter)
                            tempg = squeeze(max(green(:,3:end-3,:,volumeCounter),[],2))';
                            tempr = squeeze(max(red(:,3:end-3,:,volumeCounter),[],2))';
                            imwrite(uint16(tempg), imgToSave_topMIP_G, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                            imwrite(uint16(tempr), imgToSave_topMIP_R, 'tif', 'Compression', 'none', 'WriteMode', 'Append');

                            rgb(:,:,1) = imadjust(uint8(255*(double(-80+(tempr))./double(max(max(-80+(tempr)))))),[0.1 0.7]);
                            rgb(:,:,2) = imadjust(uint8(255*(double(-100+(tempg))./double(max(max(-100+(tempg)))))),[0.1 0.7]);
                            rgb(:,:,3) = zeros(size(tempg));
%                             rgb(:,:,2) = rgb(:,:,2)*0.8;
                            rgb = rgb*1.1;

                            ss = size(rgb);
                            tempy = interp3(conversionFactors(1)*[1:ss(2)],conversionFactors(3)*[1:ss(1)]',[1:ss(3)],double(rgb),[1:ss(2)*conversionFactors(1)],[1:ss(1)*conversionFactors(3)]',[1:ss(3)]);
                            ss = size(tempy);
                           
                            tempy = cat(1,tempy,256*ones(2,ss(2),3));
                            
                            tempgs = squeeze(max(green(:,3:end-3,:,volumeCounter),[],3))';
                            temprs = squeeze(max(red(:,3:end-3,:,volumeCounter),[],3))';

                            imwrite(uint16(tempgs), imgToSave_sideMIP_G, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                            imwrite(uint16(tempgs), imgToSave_sideMIP_R, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                            rgb2(:,:,1) = imadjust(uint8(255*(double(-80+(temprs))./double(max(max(-80+(temprs)))))),[0.1 0.9]);
                            rgb2(:,:,2) = imadjust(uint8(255*(double(-100+(tempgs))./double(max(max(-100+(tempgs)))))),[0.1 0.9]);
                            rgb2(:,:,3) = zeros(size(tempgs));
                            rgb2 = rgb2*1.1; 
                            ss = size(rgb2);

                            tempy2 = interp3(conversionFactors(1)*[1:ss(2)],conversionFactors(2)*[1:ss(1)]',[1:ss(3)],double(rgb2),[1:ss(2)*conversionFactors(1)],[1:ss(1)*conversionFactors(2)]',[1:ss(3)]);
                            
                            tempy = cat(1,tempy, tempy2);
                            %             disp(volumeCounter)
                        else%log vid
                            tempg = log(double(squeeze(max(green(:,3:end-3,:,volumeCounter),[],2))'));
                            tempr = log(double(squeeze(max(red(:,3:end-3,:,volumeCounter),[],2))'));
                            imwrite(uint16(tempg), imgToSave_topMIP_G, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                            imwrite(uint16(tempr), imgToSave_topMIP_R, 'tif', 'Compression', 'none', 'WriteMode', 'Append');

%                             rgb(:,:,1) = imadjust(uint8(255*(double(-80+(tempr))./double(max(max(-80+(tempr)))))),[0.1 0.9]);
                            rgb(:,:,1) = imadjust(uint8(255*(double(-log(90)+(tempr))./double(max(max(-log(90)+(tempr)))))),[0.1 0.7]);
                            rgb(:,:,2) = imadjust(uint8(255*(double(-log(100)+(tempg))./double(max(max(-log(100)+(tempg)))))),[0.1 0.7]);
                            rgb(:,:,3) = zeros(size(tempg));
                            rgb = rgb*1.1;

                            tempy = rgb;
                            ss = size(rgb);
                            tempy = interp3(conversionFactors(1)*[1:ss(2)],conversionFactors(3)*[1:ss(1)]',[1:ss(3)],double(rgb),[1:ss(2)*conversionFactors(1)],[1:ss(1)*conversionFactors(3)]',[1:ss(3)]);
                            ss = size(tempy);
                           
                            tempy = cat(1,tempy,256*ones(2,ss(2),3));
                            
                            
                            tempgs = log(double(squeeze(max(green(:,3:end-3,:,volumeCounter),[],3))'));
                            temprs = log(double(squeeze(max(red(:,3:end-3,:,volumeCounter),[],3))'));

                            imwrite(uint16(tempgs), imgToSave_sideMIP_G, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                            imwrite(uint16(tempgs), imgToSave_sideMIP_R, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                            rgb2(:,:,1) = imadjust(uint8(255*(double(-log(90)+(temprs))./double(max(max(-log(90)+(temprs)))))),[0.1 0.7]);
                            rgb2(:,:,2) = imadjust(uint8(255*(double(-log(100)+(tempgs))./double(max(max(-log(100)+(tempgs)))))),[0.1 0.7]);
                            rgb2(:,:,3) = zeros(size(tempgs));
                            rgb2 = rgb2*1.1;
                            ss = size(rgb2);
                            tempy2 = interp3(conversionFactors(1)*[1:ss(2)],conversionFactors(2)*[1:ss(1)]',[1:ss(3)],double(rgb2),[1:ss(2)*conversionFactors(1)],[1:ss(1)*conversionFactors(2)]',[1:ss(3)]);
                                                    
                            tempy = cat(1,tempy, tempy2);
                            %             disp(volumeCounter)
                        end
                        
                        M(volumeCounter-stt+1) = im2frame(uint8(tempy));
                    end
                    writeVideo(vidObj,M);
                    close(vidObj)
                    
                    disp('done')
                end
            end
        end
    end
    
    
    if isfield(handles,'tiff')
        if handles.tiff == 1
            dataDirectory = get(handles.foldername,'String');
            if(isdir([dataDirectory,'/tiff_stacks']) == 0)
                [dataDirectory,'/tiff_stacks']
                mkdir([dataDirectory,'/tiff_stacks']);
            end
            
            if(isdir([dataDirectory,'/tiff_stacks/',handles.scanName]) == 0)
                mkdir([dataDirectory,'/tiff_stacks/',handles.scanName]);
            end
            
            filePath = [dataDirectory,'/tiff_stacks/',handles.scanName,'/'];
            
            tic
            disp('Start Saving Scan');
             ss = size(green);
%             ss = size(SCAPE_data);
%             conversionFactors = [handles.info.info.cal.ylat; handles.info.info.cal.zdep; handles.info.info.cal.xwid];
            if numel(ss)== 3 % no for HR yet
%                 imgToSave = [filePath, 'unCorrected_' handles.scanName '.tiff'];
%                 fileCounter = 1;
%                 while (2 == exist(imgToSave, 'file'))
%                     imgToSave = [filePath, 'unCorrected_' handles.scanName '_' num2str(fileCounter) '.tiff'];
%                     fileCounter = fileCounter + 1;
%                 end
%                 
%                 if get(handles.checkbox2,'Value')
%                     
%                     % Coordinate System Correction
%                     SCAPE_data = flip(flip(SCAPE_data(:,:,:),2),3);
%                     SCAPE_data = permute(SCAPE_data, [3 1 2]);
%                     R = imref3d(size(SCAPE_data), conversionFactors(3), conversionFactors(1), conversionFactors(2));
%                     
%                     % Correct for skew
%                     
%                     %                        figure(99); subplot(1,2,1)
%                     %                        imagesc([0:size(SCAPE_data, 1)-1]*conversionFactors(3),[0:size(SCAPE_data, 3)-1]*conversionFactors(2), log(double(squeeze(max(SCAPE_data(:,round(end/3):round(2*end/3), :), [], 2))))')
%                     %                        colormap gray; axis image
%                     
%                     delta = str2num(get(handles.SkewAng,'String'))*pi/180;
%                     
%                     affineMatrix = [1 0 0 0;
%                         0 1 0 0;
%                         0 cot(delta) 1 0;
%                         0 0 0 1];
%                     tform = affine3d(affineMatrix);
%                     [SCAPE_data, ~] = imwarp(SCAPE_data, R, tform);
%                     
%                     %                     figure(99); subplot(1,2,2)
%                     %                     imagesc([0:size(SCAPE_data, 1)-1]*conversionFactors(3),[0:size(SCAPE_data, 3)-1]*conversionFactors(2), log(double(squeeze(max(SCAPE_data(:,round(end/3):round(2*end/3), :), [], 2))))')
%                     %                     colormap gray; axis image
%                     imgToSave = [filePath, handles.scanName '.tiff'];
%                     disp('Skew correction applied')
%                 end
%                 
%                 for j = 1:size(SCAPE_data, 3)
%                     imwrite(uint16(SCAPE_data(:, 1:end, j)), imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
%                 end
%                 
            else
            parfor_progress(ss(4));          
                parfor (i  =  1:ss(4),8) %scans
                    imgToSave = [filePath, 'unCorrected_' handles.scanName '_t' num2str(i) '.tiff'];
                    fileCounter = 1;
                    while (2 == exist(imgToSave, 'file'))
                        imgToSave_G = [filePath, 'Corrected_' handles.scanName '_' num2str(fileCounter) '_t' num2str(i) '_G.tiff'];
                        imgToSave_R = [filePath, 'Corrected_' handles.scanName '_' num2str(fileCounter) '_t' num2str(i) '_R.tiff'];
                        fileCounter = fileCounter + 1;
                    end
                    temp_G = green(:,:,:,i);
                    temp_R = red(:,:,:,i);
                    if get(handles.checkbox2,'Value')
                        % Coordinate System Correction
%                         temp = flip(flip(temp,2),3);
%                         temp = permute(temp, [3 1 2]);
%                         R = imref3d(size(temp), conversionFactors(3), conversionFactors(1), conversionFactors(2));
%                         
%                         % Correct for skew
%                         %                     figure(1)
%                         %                     imagesc([0:size(SCAPE_data, 1)-1]*conversionFactors(3), ...
%                         %                         [0:size(SCAPE_data, 3)-1]*conversionFactors(2), squeeze(max(SCAPE_data(:, 300:350, :), [], 2))')
%                         %                     axis image
%                         delta = str2num(get(handles.SkewAng,'String'))*pi/180;
%                         %delta = pi/2;
%                         affineMatrix = [1 0 0 0;
%                             0 1 0 0;
%                             0 cot(delta) 1 0;
%                             0 0 0 1];
%                         tform = affine3d(affineMatrix);
%                         [temp, ~] = imwarp(temp, R, tform);
                        imgToSave_G = [filePath, handles.scanName '_' num2str(i) '_G.tiff'];
                        imgToSave_R = [filePath, handles.scanName '_' num2str(i) '_R.tiff'];
                    end
                    for j = 1:size(temp_G, 3)
                        imwrite(uint16(temp_G(:, 1:end, j)), imgToSave_G, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                        imwrite(uint16(temp_R(:, 1:end, j)), imgToSave_R, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                    end
               parfor_progress;
                end
              parfor_progress(0);
                
            end
            disp('Finished Saving Scan');
        end
    end
    
    [handles] = listbox1_Callback(hObject, eventdata, handles);
    
    clear SCAPE_data;
    
else
    disp(sprintf('NOTHING TO LOAD: %s',handles.info.info.scanName));
    keyboard
end
