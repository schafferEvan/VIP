function varargout = SCAPE_simpleGalvoGUI_v13(varargin)
% SCAPE_SIMPLEGALVOGUI_V13 MATLAB code for SCAPE_simpleGalvoGUI_v13.fig
% %      SCAPE_SIMPLEGALVOGUI_V13, by itself, creates a new SCAPE_SIMPLEGALVOGUI_V13 or raises the existing
% %      singleton*.
%
%      H = SCAPE2_GALVO_GUI_V15 returns the handle to a new SCAPE_SIMPLEGALVOGUI_V13 or the handle to
%      the existing singleton*.
%
%      SCAPE_SIMPLEGALVOGUI_V13('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCAPE_SIMPLEGALVOGUI_V13.M with the given input arguments.
%
%      SCAPE_SIMPLEGALVOGUI_V13('Property','Value',...) creates a new SCAPE_SIMPLEGALVOGUI_V13 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
% %      applied to the GUI before SCAPE_simpleGalvoGUI_v13_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
% %      stop.  All inputs are passed to SCAPE_simpleGalvoGUI_v13_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDAT4 GUIHANDLES

% Edit the above text to modify the response to help SCAPE_simpleGalvoGUI_v13

% Last Modified by GUIDE v2.5 09-Mar-2018 14:03:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SCAPE_simpleGalvoGUI_v13_OpeningFcn, ...
    'gui_OutputFcn',  @SCAPE_simpleGalvoGUI_v13_OutputFcn, ...
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

% --- Executes just before SCAPE_simpleGalvoGUI_v13 is made visible.
function SCAPE_simpleGalvoGUI_v13_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SCAPE2_  galvo_GUI_v14 (see VARARGIN)

% Find path containing GUI and add it to path.
handles.GUI_filepath = mfilename('fullpath');
tmp = find(handles.GUI_filepath=='\');
handles.GUI_filepath = handles.GUI_filepath(1:tmp(end));
addpath(handles.GUI_filepath);

% insert picture in GUI
guiimg = imread([handles.GUI_filepath '/guifig2.jpg']);
axes(handles.guifig1)
image(guiimg)
axis off
axis image

axes(handles.guifig2)
image(guiimg)
axis off
axis image


%Create tab group will call scan calculations when new tab is selected by user
handles.tgroup = uitabgroup('Parent', handles.TabParent,'TabLocation', 'top','SelectionChangedFcn',@(hObject,eventdata)SCAPE_simpleGalvoGUI_v13('tgroup_SelectionChangedFcn',hObject,eventdata,guidata(hObject)) );
handles.tab1 = uitab('Parent', handles.tgroup, 'Title', 'HR');
handles.tab2 = uitab('Parent', handles.tgroup, 'Title', 'Fast');

%Place panels into each tab
set(handles.HRpanel,'Parent',handles.tab1)
set(handles.fastpanel,'Parent',handles.tab2)

%Reposition each panel to same location as panel 1
set(handles.fastpanel,'position',get(handles.HRpanel,'position'));

% create data acquisition session for dio
handles.dioSession = daq.createSession('ni');

% add digital lines to session
% 1. Trigger Line
handles.dioSession.addDigitalChannel('Dev1', 'Port0/Line0', 'OutputOnly');
% Digital shutter lines for blue (line 2) 
handles.dioSession.addDigitalChannel('Dev1', 'Port0/Line2:3', 'OutputOnly');


% STIM
% Create AO session for stim out (for those system with stim capabilities -
% usually AO0)

handles.aoStimSession = daq.createSession('ni');
handles.aoStimSession.addAnalogOutputChannel('Dev1', 'ao0', 'Voltage');
% handles.aoStimSession.addTriggerConnection('external', 'Dev1/PFI0', 'StartTrigger');
% handles.aoStimSession.addTriggerConnection('external', 'Dev1/PFI0', 'StartTrigger');

%handles.aoStimSession.Connections(1).TriggerCondition = 'RisingEdge'; %***


% Create AO session for Galvo out
% create data acquisition session for AO for scanning
handles.aoSession = daq.createSession('ni');
% create AO object for scanning adding analog output channel 0.
% Line 1 is usually Galvo in 
handles.aoSession.addAnalogOutputChannel('Dev1', 'ao1', 'Voltage');
% set AO for scanning to external triggering
handles.aoSession.addTriggerConnection('external', 'Dev1/PFI0', 'StartTrigger');
handles.aoSession.Connections(1).TriggerCondition = 'RisingEdge';
% set AO for scanning to external clock source on PFI9
handles.aoSession.addClockConnection('external', 'Dev1/PFI9', 'ScanClock');

% Create AO session for laser power control
% This is a static analog output value. No need to use
handles.aoLaserSession = daq.createSession('ni');
% Line 0 is the blue laser power + the red laser power
handles.aoLaserSession.addAnalogOutputChannel('Dev1', 'ao2', 'Voltage');


% % motorized stage - trigger for start scan (for systems with motorized
% stage capabilities)
% handles.dioStageSession = daq.createSession('ni');
% handles.dioStageSession.addDigitalChannel('Dev1', 'Port0/Line4', 'OutputOnly');
% handles.dioStageSession.addTriggerConnection('external', 'Dev1/PFI0', 'StartTrigger');
% handles.dioStageSession.Connections(1).TriggerCondition = 'RisingEdge';
% handles.dioStageSession.addClockConnection('external', 'Dev1/PFI9', 'ScanClock');

% % create data acquisition session for external exposure signal
% handles.aoExposureSession = daq.createSession('ni');
% % Add a counter channel to Channel ID 0 aka ctr0. Forces that channel to
% % output an
% handles.aoExposureSession.addCounterOutputChannel('Dev2', 0, 'PulseGeneration');
% % Set external exposure signal to external triggering
% %handles.aoSession.addTriggerConnection('external', 'Dev2/PFI0', 'StartTrigger');
% handles.aoExposureSession.Connections(1).TriggerCondition = 'RisingEdge';

% create data acquisition session for AI
handles.aiSession = daq.createSession('ni');
% create AI object adding analog input channels 0:2
% note: number of input channels is set in GUI when stimulus box is checked
aiChannels = str2num( get(handles.ai_input_channels, 'String'));
handles.aiSession.addAnalogInputChannel('Dev1', [0:aiChannels-1], 'Voltage');
% set AI to external triggering
handles.aiSession.addTriggerConnection('external', 'Dev1/PFI0', 'StartTrigger');
handles.aiSession.Connections(1).TriggerCondition = 'RisingEdge';


% set AI sample rate to default
aiSampleRate = str2double( get(handles.ai_sample_rate, 'String'));
if (aiSampleRate > handles.aiSession.RateLimit(2))
    handles.aiSession.Rate = handles.aiSession.RateLimit(2);
    errordlg(['You asked for a faster sample rate than this board can take (max = ' num2str(handles.aiSession.RateLimit(2)) ' samples per second. Setting to max.'],...
        'ERROR', modal);
elseif (aiSampleRate <= handles.aiSession.RateLimit(2))
    handles.aiSession.Rate = aiSampleRate;
end

% set AI scan length to default
handles.aiSession.DurationInSeconds = 1 + str2double( get(handles.scan_length, 'String'));

% put default galvo position at 0
handles.currentGalvoPosition = 0;

% set galvo conversion factor for input voltage to optical scan angle
handles.galvoConversionFactor = 0.25;

% Start serial communication with filter wheel 
try
    handles.filterWheelObj = serial('COM6', 'BaudRate', 115200, 'DataBits', 8, ...
        'Terminator', 'CR');
    fopen(handles.filterWheelObj);
    fprintf(handles.filterWheelObj, 'pos=1');
catch
    try
    obj = instrfind('Name', 'Serial-COM6');
    delete(obj)
    handles.filterWheelObj = serial('COM6', 'BaudRate', 115200, 'DataBits', 8, ...
        'Terminator', 'CR');
    fopen(handles.filterWheelObj);
    fprintf(handles.filterWheelObj, 'pos=1');
    catch
         h1 = msgbox('Unable to establish connection with filter wheel. Check if device is turned on and try again.');
    end
end


% Get calibration factors
numObjectives = size(get(handles.objective_list, 'String'),1);
handles.calibration_filepath = [handles.GUI_filepath 'SCAPE_calibrations.txt'];
if exist(handles.calibration_filepath,'file')==2
    fid = fopen(handles.calibration_filepath,'r');
    a = fgetl(fid); % Skip 1st line
    a = fgetl(fid); % Skip 2nd line
    
    handles.savedConversionFactors = struct([]);
    for i = 1 %:numObjectives
        calFactors = fgetl(fid);
        p = strfind(calFactors,'	');
        handles.savedConversionFactors(i).y_umPerPix = str2num(calFactors(p(2)+1:p(3)-1));
        handles.savedConversionFactors(i).z_umPerPix = str2num(calFactors(p(3)+1:p(4)-1));
        handles.savedConversionFactors(i).xK_umPerVolt = str2num(calFactors(p(4)+1:end));
    end
fclose(fid);
else
     errordlg('Cannot find SCAPE Calibration file');
     handles.savedConversionFactors(1).y_umPerPix = 1.39;
     handles.savedConversionFactors(1).z_umPerPix = 1.11;
     handles.savedConversionFactors(1).xK_umPerVolt = 313.51; 
end


% I change this parameter when I initially set up the system. (um or volts)
handles.globalGalvoOffset = 0;

% Set up listener for when the frame rate update config file changes
handles.frameRateFile = System.IO.FileSystemWatcher(handles.GUI_filepath);
handles.frameRateFile.Filter = 'cameraConfig_frameRateUpdate.txt';
handles.frameRateFile.EnableRaisingEvents = true;

% Set initial parameter variables
handles.cameraFramerate = str2num(get(handles.camera_framerate, 'String'));
handles.cameraFramerate2 = str2num(get(handles.camera_framerate2, 'String'));
handles.originalFPS = 10;
handles.originalFPS2 = 1370;
handles.framerate_newInput = 1;

%HR scan parameters
handles.scanRate = str2num(get(handles.volumetric_scan_rate, 'String'));
handles.scanFOV = str2num(get(handles.scan_FOV, 'String'));
handles.numScanSteps = str2num(get(handles.num_scan_steps, 'String'));
handles.stepSize = str2num(get(handles.step_size, 'String'));
handles.scanLength = str2num(get(handles.scan_length, 'String'));

%fast scan parameters
handles.scanRate2 = str2num(get(handles.volumetric_scan_rate2, 'String'));
handles.scanFOV2 = str2num(get(handles.scan_FOV2, 'String'));
handles.numScanSteps2 = str2num(get(handles.num_scan_steps2, 'String'));
handles.stepSize2 = str2num(get(handles.step_size2, 'String'));
handles.scanLength2 = str2num(get(handles.scan_length2, 'String'));
handles.info.loadParameters.maxFrameRate = get(handles.max_frame_rate, 'Value');

% original filter wheel value set to 'open'
handles.originalfiltval = 1;
% set alignment mode filter wheel setting to 1%
set(handles.alignpower, 'Value', 10)

% default to calculating a new framerate + updating parameters
handles.framerate_newFPS = 0;
handles.paramset = 1;

% Choose default command line output for SCAPE_simpleGalvoGUI_v13
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% find data directory
browse_for_data_directory_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SCAPE_simpleGalvoGUI_v13 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = SCAPE_simpleGalvoGUI_v13_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function scan_calculations(hObject, handles)
disp(' ');
disp('Updating...')
% Set data+scan directory
info = handles.info;
info.dataDirectory = get(handles.data_directory, 'String');
scanName_prefix = get(handles.scan_name_stem, 'String');

infoFileList = dir([info.dataDirectory, '*.mat']);
fileCounter = 1;
for i = 1:length(infoFileList)
    tmp = infoFileList(i).name;
    if (~isempty(strfind(tmp, [scanName_prefix '_run'])))
        fileCounter = fileCounter+1;
    end
end

info.scanName = [scanName_prefix '_run' num2str(fileCounter)];

%adds suffix
if(~isempty(get(handles.scan_name_suffix, 'String')))
    info.scanName = [info.scanName '_' get(handles.scan_name_suffix, 'String')];
end
%adds HR tag
%kpedit - check if HR tab is selected
%isHR = get(handles.HR_checkbox, 'Value');

isHR = strcmp('HR',handles.tgroup.SelectedTab.Title);
if(isHR == 1)
    info.scanName = [info.scanName '_HR'];
end

% Step 1: Set ROI values on GUI
isFullFrame = get(handles.fullFrame, 'Value');
if (isFullFrame == 0)
    if(isHR == 1)
        info.camera.xROI = str2num(get(handles.x_ROI, 'String'));
        info.camera.yROI = str2num(get(handles.y_ROI, 'String'));
        info.camera.x_left = str2num(get(handles.camera_left, 'String'));
        temp1 = get(handles.pixel_binning, 'String');
        stimPattern = get(handles.pixel_binning, 'Value');
        temp3 = cell2mat(temp1(stimPattern));
        info.camera.binFactor = str2num(temp3(1));
        clear temp1 temp2 temp3
    else
        info.camera.xROI = str2num(get(handles.x_ROI2, 'String'));
        info.camera.yROI = str2num(get(handles.y_ROI2, 'String'));
        info.camera.x_left = str2num(get(handles.camera_left2, 'String'));
        temp1 = get(handles.pixel_binning2, 'String');
        stimPattern = get(handles.pixel_binning2, 'Value');
        temp3 = cell2mat(temp1(stimPattern));
        info.camera.binFactor = str2num(temp3(1));
        clear temp1 temp2 temp3
    end
else  
    info.camera.xROI = 2042;
    info.camera.yROI = 2042;
    info.camera.x_left = 1;
    if(isHR == 1)
        temp1 = get(handles.pixel_binning, 'String');
        stimPattern = get(handles.pixel_binning, 'Value');
    else
        temp1 = get(handles.pixel_binning2, 'String');
        stimPattern = get(handles.pixel_binning2, 'Value');
    end
    temp3 = cell2mat(temp1(stimPattern));
    info.camera.binFactor = str2num(temp3(1));
    clear temp1 temp2 temp3
    handles.paramset = 1; % will update paramters on Solis
end

% Step 2: Check objective used and set calibration factors
X = get(handles.objective_list, 'String');
value = get(handles.objective_list, 'Value');

if (iscell(X))
    handles.objectiveUsed = cell2mat(X(value));
else
    handles.objectiveUsed = X;
end
info.objective = handles.objectiveUsed;
clear X value;
if (isempty(strfind(handles.objectiveUsed, '10x')) == 0)
    cFact = handles.savedConversionFactors(1);
    info.GUIcalFactors.z_umPerPix = cFact.z_umPerPix;     % Does change as a function of objective magnification
    info.GUIcalFactors.y_umPerPix = cFact.y_umPerPix;   % Does change as a function of objective magnification
    info.GUIcalFactors.xK_umPerVolt = cFact.xK_umPerVolt; % Doesn't change as a function of objective magnification (Volts to galvo/um)
elseif(isempty(strfind(handles.objectiveUsed, '20x')) == 0)
    if (isempty(strfind(handles.objectiveUsed, 'Olympus'))==0)
        cFact = handles.savedConversionFactors(1);
        scaleFactor = 7.77777777777777/3.5;
        info.GUIcalFactors.z_umPerPix = cFact.z_umPerPix/scaleFactor;   % Does change as a function of objective magnification
        info.GUIcalFactors.y_umPerPix = cFact.y_umPerPix/scaleFactor;   % Does change as a function of objective magnification
        info.GUIcalFactors.xK_umPerVolt = cFact.xK_umPerVolt; % Doesn't change as a function of objective magnification (Volts to galvo/um)
        
    elseif (isempty(strfind(handles.objectiveUsed, 'Nikon'))==0)
        cFact = handles.savedConversionFactors(1);
        scaleFactor = 7/3.5;
        info.GUIcalFactors.z_umPerPix = cFact.z_umPerPix/scaleFactor;     % Does change as a function of objective magnification
        info.GUIcalFactors.y_umPerPix = cFact.y_umPerPix/scaleFactor;   % Does change as a function of objective magnification
        info.GUIcalFactors.xK_umPerVolt = cFact.xK_umPerVolt; % Doesn't change as a function of objective magnification (Volts to galvo/um)
    end
end

% Step 3: Set initial scan parameters
% Camera Frame Rate
inPreview = get(handles.preview, 'Value');
if(isHR == 1)
    handles.cameraFramerate = str2num(get(handles.camera_framerate, 'String'));
    info.camera.framerate = handles.cameraFramerate;
else
    handles.cameraFramerate2 = str2num(get(handles.camera_framerate2, 'String'));
    info.camera.framerate = handles.cameraFramerate2;
end
if (isFullFrame == 0 && inPreview == 0)
    if (handles.framerate_newInput == 1)
        disp('New fps input!!!')
        if(isHR == 1)
            info.loadParameters.desiredCameraFramerate = handles.cameraFramerate;
        else
            info.loadParameters.desiredCameraFramerate = handles.cameraFramerate2;
        end
        handles.framerate_newInput = 0;
    end
else
    disp('fps input set to 10!')
    info.loadParameters.desiredCameraFramerate = 10;
end

if(isHR == 1)
    % Scan Rate (Volumes/sec)
    info.daq.scanRate = handles.scanRate;
    % Scan FOV:
    info.daq.scanFOV = handles.scanFOV;
else
    % info.camera.framerate = str2num(get(handles.camera_framerate, 'String'));
    info.loadParameters.maxFrameRate = get(handles.max_frame_rate, 'Value');
    % Scan Rate (Volumes/sec)
    info.daq.scanRate = handles.scanRate2;
    % Scan FOV:
    info.daq.scanFOV = handles.scanFOV2;
end
offset_um  = str2num(get(handles.galvo_um_offset, 'String'));


% Number of frames allocated to each volume (add an extra frame for the
% flyback)

if(isHR == 1)
    info.daq.pixelsPerLine = handles.numScanSteps+1;
else
    info.daq.pixelsPerLine = handles.numScanSteps2+1;
end


% Step 4: Set information regarding scan duration
% Total number of scans to acquire and number of seconds to acquire

if (isHR == 1)
    info.daq.numberOfScans = 1;
    info.daq.scanLength = (info.daq.pixelsPerLine/info.camera.framerate);
    info.HR = 1;
    %   info.daq.numberOfScans = 1;
    set(handles.scan_length, 'Enable', 'off');
    set(handles.scan_length, 'BackgroundColor', [0.5 0.5 0.5]);
    set(handles.scan_length, 'String', num2str(info.daq.scanLength));
else
    timePerVolume = (info.daq.pixelsPerLine/info.camera.framerate);
    scanLength = handles.scanLength2;
    info.daq.numberOfScans = floor(scanLength/timePerVolume);
    if ((scanLength/timePerVolume)<1)
        if (inPreview==0)    
        disp('Must take at least 1 volume');
        info.daq.numberOfScans = 1;
        end
    end
    handles.scanLength2  =  info.daq.numberOfScans*timePerVolume;
    info.daq.scanLength = handles.scanLength2;
    set(handles.scan_length2, 'Enable', 'on');
    set(handles.scan_length2, 'BackgroundColor', [1 1 1]);
    set(handles.scan_length2, 'String', num2str(info.daq.scanLength));
end

% Number of frames acquired by camera
info.camera.kineticSeriesLength =  num2str(ceil(info.daq.scanLength*info.camera.framerate)+100);


% Step 5: Set up scan waveform
info.daq.sawtooth = 1;
info.daq.scanAngle = (info.daq.scanFOV/info.GUIcalFactors.xK_umPerVolt);
info.daq.galvoOffset = offset_um/info.GUIcalFactors.xK_umPerVolt;
globalOffset_volts = handles.globalGalvoOffset/info.GUIcalFactors.xK_umPerVolt;
scanPattern = zeros(info.daq.pixelsPerLine, 1);
maxVal = globalOffset_volts+info.daq.galvoOffset+info.daq.scanAngle/2;
minVal = globalOffset_volts+info.daq.galvoOffset-info.daq.scanAngle/2;
scanPattern(1:end-1) = linspace(minVal, maxVal, info.daq.pixelsPerLine-1);
scanPattern(end) = minVal;
scanPattern = repmat(scanPattern, [info.daq.numberOfScans, 1]);
info.GUIcalFactors.x_umPerPix = info.GUIcalFactors.xK_umPerVolt*info.daq.scanAngle/(info.daq.pixelsPerLine-2);

% Step 6: Set up stim waveform
info.loadParameters.useStim = get(handles.use_stim_checkbox, 'Value');
t = ([1:length(scanPattern)]-1)/info.camera.framerate;
diffTime = mean(diff(t));

if ((info.loadParameters.useStim==1)&&(strcmp('HR', handles.tgroup.SelectedTab.Title)==0))
    stimPattern     = scanPattern*0;
    stimParameters  = get(handles.stimTable, 'data');
    info.daq.stim.Vout  = 0;
    info.daq.stim.parameters = stimParameters;
    
    stimOnsets = stimParameters(:, 1);
    tmp = strcmp(stimOnsets, '');
    stimDurations = stimParameters(:, 2);
    tmp = tmp+strcmp(stimDurations, '');
    stimVoltage = stimParameters(:, 3);
    tmp = tmp+strcmp(stimVoltage, '');
    numStims = sum(uint16(tmp==0));
    for i = 1:numStims
        onset_sec = str2num(stimOnsets{i});
        duration_sec = str2num(stimDurations{i})/1000;
        tmp = abs(t-onset_sec);
        ind = find(tmp == min(tmp));
        indDuration = round(duration_sec/diffTime);
        voltage = str2num(stimVoltage{i});
        stimPattern(ind:ind+indDuration-1) = voltage;
    end
    %
    stimPattern = stimPattern(1:length(scanPattern));
    plot(handles.stim_axes, t, stimPattern);
    xlim([t(1) t(end)]);
else
    stimPattern = scanPattern*0;
    info.daq.stim.Vout = 0;
    cla(handles.stim_axes);
end

% STIM - Stim out
handles.info.daq.stim.Vout = stimPattern';
info.daq.scanWaveform = scanPattern';


% Step 7: Read laser parameters
useBlueLaser = get(handles.blue_laser_checkbox, 'Value');
useRedLaser = get(handles.red_laser_checkbox, 'Value');
info.shutterToggle = [useBlueLaser useRedLaser];
%Modifying laser power info to reflect preview mode changes
if(inPreview == 0)
    info.blue_laser_output_power = str2num(get(handles.blue_laser_output_power, 'String'));
    info.red_laser_output_power = str2num(get(handles.red_laser_output_power, 'String'));
else
    if(str2num(get(handles.blue_laser_output_power, 'String'))>2)
        info.blue_laser_output_power = 2;
    else
        info.blue_laser_output_power = 0;
    end
    if(str2num(get(handles.red_laser_output_power, 'String'))>2)
        info.red_laser_output_power = 2;
    else
        info.red_laser_output_power = 0;
    end
    
end

% Step 8: Update AI channels, sample rate and scan duration
info.daq.aiChannels = [0:str2num( get(handles.ai_input_channels, 'String'))-1];

%this throws errors sometimes...
try
while 1 <= length(handles.aiSession.Channels)
    handles.aiSession.removeChannel(1);
end
handles.aiSession.addAnalogInputChannel('Dev1', info.daq.aiChannels, 'Voltage');
info.daq.aiSampleRate = str2double( get(handles.ai_sample_rate, 'String'));
handles.aiSession.Rate = info.daq.aiSampleRate;
handles.aiSession.DurationInSeconds = 1 + info.daq.scanLength;
catch
    disp('AI not updated')
end

% Step 9:
set(handles.data_save_path, 'String', [info.dataDirectory info.scanName]);

% Step 10: Update the laser power
filterstring = get(handles.filter_wheel, 'String');
laserOD = get(handles.filter_wheel, 'Value');
info.laser_power = cell2mat(filterstring(laserOD));

% Step 11: Update experiment notes
info.experiment_notes = get(handles.experiment_notes, 'String');

% Step 12: update handles
handles.info = info;
handles.framerateChanged = addlistener(handles.frameRateFile,'Changed', ...
    @(hObject,eventdata)SCAPE_simpleGalvoGUI_v13('frame_rate_changed',hObject,eventdata,handles));

startpreview = get(handles.solispreview, 'Value');
% Step 13: Write the config file for Zyla
writeConfigFile(hObject, handles);

% update solis only if preview is on and parameters have changed OR framerate needs to be recalculated
if (startpreview && handles.paramset == 1) || handles.framerate_newFPS
    system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_stopPreview.ahk');
    system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_setParameters.ahk');
    pause(2); % Needs this to give set parameters time to run (probably to update listener as well);
    if startpreview == 1
        system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_startPreview.ahk');
    end
    system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\returnToGUI.ahk');
    handles.paramset = 0; handles.framerate_newFPS=0;  % reset b/c no new parameters need updating in solis
end

guidata(hObject, handles);

disp('Scan Calculations Complete');

function scan_Callback(hObject, eventdata, handles)
% get motorized stage ready for movement
isHR = strcmp('HR', handles.tgroup.SelectedTab.Title);
if (get(handles.motorstage, 'Value') && isHR==0)
    pos = handles.stage.GetPosition_Position(0);
    newpos = str2num(get(handles.stage_dist, 'String'))+ pos;
    if newpos<0; newpos = 0; end
    if newpos>50; newpos = 50; end
    handles.stage.SetAbsMovePos(0, single(newpos));
    handles.stage.SetKCubeTriggerParams(0, 3, 1, 12, 1) %will move to abs position on trigger in
end

handles.paramset = 0; %run scan calc without updating solis yet
scan_calculations(hObject,handles);

system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_stopPreview.ahk'); %precautionary step
system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_setParameters.ahk');
pause(2);
system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_startAcquisition.ahk');
pause(1.5)
system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\returnToGUI.ahk');
clc
cla(handles.stim_axes);

handles.paramset = 0; handles.framerate_newFPS = 0; %will not update solis
scan_calculations(hObject,handles);
handles = guidata(hObject);

% reset handles.stopButtonFlag
handles.stopButtonFlag = 0;

% move galvo to start pos
handles.aoSession.outputSingleScan(handles.info.daq.scanWaveform(1)')


% STIM - set stim line to 0
if (1 == get(handles.use_stim_checkbox, 'Value'))
    handles.aoStimSession.outputSingleScan(handles.info.daq.stim.Vout(1));
end

% set digital line value to LOW keeping shutters closed
handles.dioSession.outputSingleScan([0 0 0]);

% % motorized stage
% if (get(handles.motorstage, 'Value') && isHR==0)
% handles.dioStageSession.outputSingleScan(0);
% end

% set scan rate
handles.aoSession.Rate = handles.info.camera.framerate;

% % STIM - set stim output rate
% handles.aoStimSession.Rate = handles.info.daq.stim.fout;

%put data on AO object for scanning
handles.aoSession.queueOutputData(handles.info.daq.scanWaveform');

% %motorized stage
% if (get(handles.motorstage, 'Value') && isHR==0)
% handles.dioStageSession.queueOutputData(1);
% end

% STIM - put data on AO object for stimming
if (1 == get(handles.use_stim_checkbox, 'Value'))
    handles.aoStimSession.queueOutputData(handles.info.daq.stim.Vout');
end

% save info file
info = handles.info;
info.scanStartTimeApprox = datestr(now);

info.camera.kineticSeriesLength = num2str(info.camera.framerate*info.daq.scanLength+100);
isHR = strcmp('HR', handles.tgroup.SelectedTab.Title);
if(isHR == 1)
    info.camera.framerate = str2num((get(handles.camera_framerate, 'String')));
else
    info.camera.framerate = str2num((get(handles.camera_framerate2, 'String')));
end
info.loadParameters.maxFrameRate = get(handles.max_frame_rate, 'Value');
info.daq.pixelrate = info.camera.framerate;
info.loadParameters.isHR = isHR;

% update AI logging file information
info.daq.logFilePath = [info.dataDirectory info.scanName '_stim_data.bin'];
info.scanStatus = 'incomplete';

% save preliminary info file
infoFilePath = [info.dataDirectory info.scanName '_info.mat'];
save(infoFilePath, 'info');
set(handles.data_save_path, 'String', [info.dataDirectory info.scanName]);

% update handles
handles.info = info;
guidata(hObject, handles);

% open log file
FID = fopen(info.daq.logFilePath, 'w', 'n');

% add listener to AI
lh = handles.aiSession.addlistener('DataAvailable', @(src, event)logData(src, event, FID));

% start AO for scan
handles.aoSession.startBackground;

% %motorized stage
% if (get(handles.motorstage, 'Value') && isHR==0)
% handles.dioStageSession.startBackground;
% end

% % start AO for stimming
% if (1 == get(handles.use_stim_checkbox, 'Value'))
%     handles.aoStimSession.startBackground;
% end

% start AI
handles.aiSession.startBackground;

% open shutters
% CHECK THIS
handles.dioSession.outputSingleScan([0 info.shutterToggle]);
handles.info.shutterToggle

% send start trigger by setting digital line to HIGH
tic;

handles.dioSession.outputSingleScan([1 info.shutterToggle]);

% wait until AO done
handles.aoSession.wait;

scanTime = toc;

% close shutters
handles.dioSession.outputSingleScan([0 0 0]);

% Unclick the buttons
set(handles.blue_shutter, 'Value', 0);
set(handles.red_shutter, 'Value', 0);
% set(handles.green_shutter, 'Value', 0)


% wait until AO stim session done
if (1 == get(handles.use_stim_checkbox, 'Value'))
        handles.aoStimSession.wait;
end

disp(['Done scan! Took: ' num2str(scanTime) ' seconds']);

% wait until AO done
handles.aoSession.wait;

% stop AI
handles.aiSession.stop;
delete(lh);
fclose(FID);


% %motorized stage
% if (get(handles.motorstage, 'Value') && isHR==0)
% handles.dioStageSession.stop;
% end


% save final info file
if (handles.aiSession.ScansAcquired ~= handles.aiSession.NumberOfScans)
    info.scanStatus = 'scan stopped by user. scan did not finish';
elseif (handles.aiSession.ScansAcquired == handles.aiSession.NumberOfScans)
    info.scanStatus = 'scan complete!';
end

% save final info file
save(infoFilePath, 'info');
set(handles.data_save_path, 'String', infoFilePath);

% update current position display
set(handles.galvo_position_display, 'String', num2str(handles.info.daq.scanWaveform(end) / handles.galvoConversionFactor));
set(handles.galvo_position_slider, 'Value', handles.info.daq.scanWaveform(end) / handles.galvoConversionFactor);

% Center the galvo and turn laser power to zero
% KPedit - added in blue laser
laser_power_r = get(handles.red_laser_output_power, 'String');
laser_power_b = get(handles.blue_laser_output_power, 'String');
set(handles.red_laser_output_power, 'String', '0');
set(handles.blue_laser_output_power, 'String', '0');
center_galvo_Callback(hObject, eventdata, handles)
set(handles.red_laser_output_power, 'String', laser_power_r);
set(handles.blue_laser_output_power, 'String', laser_power_b);

% % restore normal function to motorized stage and update positions
% if (get(handles.motorstage, 'Value') && isHR==0)
%     handles.stage.SetKCubeTriggerParams(0, 0, 1, 0, 1); %will move to abs position without trigger
%     handles.stage.EnableHWChannel(0);
%     handles.stage.StopImmediate(0); %stop motion if not already stopped
%     % Verify stage position
%     pause(1)
%     pos = handles.stage.GetPosition_Position(0);
%     set(handles.stage_slider, 'Value', pos);
%     set(handles.stage_pos, 'String', num2str(pos));  
% end


% load in stim data
FID = fopen(info.daq.logFilePath, 'r');
numChannels = length(info.daq.aiChannels);
stimData = fread(FID, [numChannels + 1, inf], 'double');
fclose(FID);
handles.lastStimData.time = stimData(1, :);
handles.lastStimData.channels = stimData(2:end, :);
guidata(hObject, handles);

% plot stim data
channelsToPlot = 1 + str2num( get(handles.ai_channels_to_plot, 'String'));
channelCount = 1;

% kpedit - do not always want to plot this
figure(101);

set(101, 'Position', [5 502 1241 420]);
for channelIdx = channelsToPlot
    subplot(1, length(channelsToPlot), channelCount);
    plot(handles.lastStimData.time, handles.lastStimData.channels(channelsToPlot(channelCount), :));
    %plot(handles.lastStimData.time(timeToPlot), handles.lastStimData.channels(channelsToPlot(channelIdx), timeToPlot));
    %plot(handles.lastStimData.channels(channelsToPlot(channelIdx), timeToPlot));
    xlim([min(handles.lastStimData.time) max(handles.lastStimData.time)]);
    
    switch channelIdx
       case 1
            title('Channel 0: Galvo Drive', 'FontSize', 12);
        case 2
            title('Channel 1: Galvo Feedback', 'FontSize', 12);
        case 3
            title('Channel 2: Frame ds', 'FontSize', 12);
        case 4
            title('Channel 3: Stimulus', 'FontSize', 12);
            ylim([-5 5])
    end
    channelCount = channelCount+1;
end

handles = rmfield(handles,'lastStimData'); % Beth added 12/18/14

scan_calculations(hObject, handles);
handles = guidata(hObject);

guidata(hObject, handles);


% Tab change - update parameters when HR or Fast tab is selected
function tgroup_SelectionChangedFcn(hObject, eventdata, handles)
handles.paramset = 1;
guidata(hObject, handles)
scan_calculations(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code + scan calculations callback %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fullFrame_Callback(hObject, eventdata, handles)
fullFrame = get(handles.fullFrame, 'Value');
inPreview = get(handles.preview, 'Value');
%kpedit - check if HR tab is selected
isHR = strcmp('HR', handles.tgroup.SelectedTab.Title);
if (fullFrame == 1)
    if (get(handles.solispreview, 'Value') == 0)
        set(handles.solispreview, 'Value', 1)
    else 
        system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_stopPreview.ahk');
        system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\returnToGUI.ahk');
    end
    %Disable preview mode
    if(inPreview == 1)
        set(handles.preview, 'Value', 0);
        %Enable laser power fields
        set(handles.blue_laser_output_power, 'Enable', 'on');
        set(handles.red_laser_output_power, 'Enable', 'on');
        %Enable framerate field and max framerate box.
        if(isHR == 1)
            set(handles.camera_framerate, 'Enable', 'on');
        else
            set(handles.camera_framerate2, 'Enable', 'on');
            set(handles.max_frame_rate, 'Enable', 'on');
        end
        
        %turn off lasers if on
        blueShutterValue = get(handles.blue_shutter, 'Value');
        redShutterValue = get(handles.red_shutter, 'Value');
        if (blueShutterValue == 1)
            set(handles.blue_shutter, 'Value', 0)
        end
        if (redShutterValue == 1)
            set(handles.red_shutter, 'Value', 0)
        end
        digitalLines = [0 0 0]; 
        handles.dioSession.outputSingleScan(digitalLines);

        %reset filter wheel
        set(handles.filter_wheel, 'Value', handles.originalfiltval);
        input = ['pos=' num2str(handles.originalfiltval)];
        fprintf(handles.filterWheelObj, input);
        previewwason = 1;
    else
        previewwason = 0;
    end
    if(isHR==1)
        set(handles.x_ROI, 'Enable', 'off');
        set(handles.y_ROI, 'Enable', 'off');
        set(handles.camera_left, 'Enable', 'off');
        %Disable framerate field
        set(handles.camera_framerate, 'Enable', 'off');
        if ~previewwason
        handles.originalFPS = handles.cameraFramerate;
        end
    else
        set(handles.x_ROI2, 'Enable', 'off');
        set(handles.y_ROI2, 'Enable', 'off');
        set(handles.camera_left2, 'Enable', 'off');
        %Disable framerate field
        set(handles.camera_framerate2, 'Enable', 'off');
        %Disable max framerate box
        set(handles.max_frame_rate, 'Enable', 'off');
        if ~previewwason
        handles.originalFPS2 = handles.cameraFramerate2;
        end
    end
    
    %Marking framerate to be recalculatedm
    handles.framerate_newInput = 1;
else
     if (get(handles.solispreview, 'Value') == 1)
     set(handles.solispreview, 'Value', 0)
     end
     if(isHR==1)
        set(handles.x_ROI, 'Enable', 'on');
        set(handles.x_ROI, 'BackgroundColor', [1 1 1]);
        set(handles.y_ROI, 'Enable', 'on');
        set(handles.y_ROI, 'BackgroundColor', [1 1 1]);
        set(handles.camera_left, 'Enable', 'on');
        set(handles.camera_left, 'BackgroundColor', [1 1 1]);
        %Enabling framerate field, max framerate box
        set(handles.camera_framerate, 'Enable', 'on');
        set(handles.camera_framerate, 'BackgroundColor', [1 1 1]);
        handles.cameraFramerate = handles.originalFPS;
        set(handles.camera_framerate, 'String', num2str(handles.originalFPS));
    else
        set(handles.x_ROI2, 'Enable', 'on');
        set(handles.x_ROI2, 'BackgroundColor', [1 1 1]);
        set(handles.y_ROI2, 'Enable', 'on');
        set(handles.y_ROI2, 'BackgroundColor', [1 1 1]);
        set(handles.camera_left2, 'Enable', 'on');
        set(handles.camera_left2, 'BackgroundColor', [1 1 1]);
        %Enabling framerate field, max framerate box
        set(handles.camera_framerate2, 'Enable', 'on');
        set(handles.camera_framerate2, 'BackgroundColor', [1 1 1]);
        set(handles.max_frame_rate, 'Enable', 'on');
        set(handles.max_frame_rate, 'BackgroundColor', [1 1 1]);
        handles.cameraFramerate2 = handles.originalFPS2;
        set(handles.camera_framerate2, 'String', num2str(handles.originalFPS2));
    end
    %Marking framerate to be recalculated
    handles.framerate_newInput = 1;
    
end

blue_laser_output_power_Callback(hObject, eventdata, handles);
red_laser_output_power_Callback(hObject, eventdata, handles);

handles.paramset = 1; %updates parameters
solispreview_Callback(hObject, eventdata, handles)


function load_info_file_Callback(hObject, eventdata, handles)
presentDirectory = pwd;
cd('E:\');
[fileName, pathName, filterIndex] = uigetfile('*.mat');
cd (pathName)
load(fileName)
cd(presentDirectory)

if info.loadParameters.isHR == 1
    handles.tgroup.SelectedTab = handles.tab1;
    % Camera Properties:
    binLookup = [1 2 3 4 8];
    tmp = find(binLookup==info.camera.binFactor);
    set(handles.pixel_binning, 'Value', tmp);
    set(handles.x_ROI, 'String', info.camera.xROI);
    set(handles.y_ROI, 'String', info.camera.yROI);
    set(handles.camera_left, 'String', info.camera.x_left)
    
    % Set frame rate
    handles.cameraFramerate = info.loadParameters.desiredCameraFramerate;
    handles.framerate_newInput = 1;
    set(handles.camera_framerate, 'String', num2str(handles.cameraFramerate));
    
    % Set volume rate
    set(handles.volumetric_scan_rate, 'String', info.daq.scanRate);
    handles.scanRate = info.daq.scanRate;
    
    % Set Run Length:
    numVolumes = info.daq.numberOfScans;
    timePerVolume = 1/handles.scanRate;
    handles.scanLength = numVolumes/handles.scanRate +(timePerVolume*0.1);
    set(handles.scan_length, 'String', num2str(handles.scanLength));
    
    % Set number of steps
    handles.numScanSteps = info.daq.pixelsPerLine-1;
    set(handles.num_scan_steps, 'String', num2str(handles.numScanSteps));
    % Set scan FOV
    handles.scanFOV = info.daq.scanFOV;
    set(handles.scan_FOV, 'String', num2str(handles.scanFOV));
    % Set galvo offset
    offset_um = info.daq.galvoOffset*info.GUIcalFactors.xK_umPerVolt;
    set(handles.galvo_um_offset, 'String', num2str(offset_um));
    scan_FOV_Callback(hObject, eventdata, handles);
    
else
    handles.tgroup.SelectedTab = handles.tab2;
    % Camera Properties:
    binLookup = [1 2 3 4 8];
    tmp = find(binLookup==info.camera.binFactor);
    set(handles.pixel_binning2, 'Value', tmp);
    set(handles.x_ROI2, 'String', info.camera.xROI);
    set(handles.y_ROI2, 'String', info.camera.yROI);
    set(handles.camera_left2, 'String', info.camera.x_left)
    
    % Set frame rate
    handles.cameraFramerate2 = info.loadParameters.desiredCameraFramerate;
    handles.framerate_newInput = 1;
    set(handles.camera_framerate2, 'String', num2str(handles.cameraFramerate2));
    set(handles.max_frame_rate, 'Value', info.loadParameters.maxFrameRate);
    
    % Set volume rate
    set(handles.volumetric_scan_rate2, 'String', info.daq.scanRate);
    handles.scanRate2 = info.daq.scanRate;
    
    % Set Run Length:
    numVolumes = info.daq.numberOfScans;
    timePerVolume = 1/handles.scanRate2;
    handles.scanLength2 = numVolumes/handles.scanRate2 +(timePerVolume*0.1);
    set(handles.scan_length2, 'String', num2str(handles.scanLength2));
    
    % Set number of steps
    handles.numScanSteps2 = info.daq.pixelsPerLine-1;
    set(handles.num_scan_steps2, 'String', num2str(handles.numScanSteps2));
    % Set scan FOV
    handles.scanFOV2 = info.daq.scanFOV;
    set(handles.scan_FOV2, 'String', num2str(handles.scanFOV2));
    % Set galvo offset
    offset_um = info.daq.galvoOffset*info.GUIcalFactors.xK_umPerVolt;
    set(handles.galvo_um_offset, 'String', num2str(offset_um));
    scan_FOV2_Callback(hObject, eventdata, handles);
    
end


% Laser power
set(handles.red_laser_output_power, 'String', num2str(info.red_laser_output_power));
set(handles.blue_laser_output_power, 'String', num2str(info.blue_laser_output_power));


% % Filter wheel 
temp = get(handles.filter_wheel, 'String');
set(handles.filter_wheel, 'Value', find(strcmp(temp, info.laser_power) == 1));
userChoice = num2str(get(handles.filter_wheel, 'Value'));
input = ['pos=' userChoice];
fprintf(handles.filterWheelObj, input);

% Set objective
temp = get(handles.objective_list, 'String');
set(handles.objective_list, 'Value', find(strcmp(temp, info.objective) == 1));

% Set stim
set(handles.use_stim_checkbox, 'Value', info.loadParameters.useStim)
if (info.loadParameters.useStim == 1)
    set(handles.stimTable, 'data', info.daq.stim.parameters);
end
scan_calculations(hObject, handles);

function filter_wheel_Callback(hObject, eventdata, handles)

userChoice = num2str(get(handles.filter_wheel, 'Value'));
input = ['pos=' userChoice];

try
    fprintf(handles.filterWheelObj, input);
catch
    clc
    try
    obj = instrfind('Name', 'Serial-COM6');
    delete(obj)
    handles.filterWheelObj = serial('COM6', 'BaudRate', 115200, 'DataBits', 8, ...
        'Terminator', 'CR');
    fopen(handles.filterWheelObj);
    catch
        h1 = msgbox('Unable to establish connection with filter wheel. Check if device is turned on and try again.');
    end
    fprintf(handles.filterWheelObj, input);
    
end

scan_calculations(hObject, handles);    


function browse_for_data_directory_Callback(hObject, eventdata, handles)
handles.info.directory = [uigetdir('E:\','Select Directory to save data to') '\'];
cd(handles.info.directory);
% For some reason, any directory chosen that is not a drive (C:, D:, etc)
% is not returned with a backslash, but if a drive is chosen, a backslash
% is appended.
if (handles.info.directory(end-1:end) == '\\')
    handles.info.directory = handles.info.directory(1:end-1);
end

if 0 == handles.info.directory
    errordlg('You did not select a directory. Please choose again.','okay')
    return;
end

set(handles.data_directory, 'String', handles.info.directory);
guidata(hObject,handles);

scan_calculations(hObject, handles);

function blue_laser_output_power_Callback(hObject, eventdata, handles)
inPreview = get(handles.preview, 'Value');

%In preview mode, laser power is set to 2 if not set to 1 or 0 prior
if (inPreview == 0)
    laserPowerBlue = str2num(get(handles.blue_laser_output_power, 'String'));
    laserPowerRed = str2num(get(handles.red_laser_output_power, 'String'));
else
    if(str2num(get(handles.blue_laser_output_power, 'String'))>=2)
        laserPowerBlue = 2;
    else
        laserPowerBlue = 0;
    end
end

% Old formula (before changing to 2000Ohm input impedance)
% voltPower = laserPower*0.0694-0.1140;
% NEW FORMULA:
voltPowerBlue =(laserPowerBlue*laserPowerBlue)*0.0001+(laserPowerBlue*0.0305)+0.0012;  %SCAPE4 lasers
%voltPowerBlue =(laserPowerBlue*laserPowerBlue)*1.423e-05+(laserPowerBlue*0.03241)+0.002877; %Uptown SCAPE lasers
voltPowerBlue(voltPowerBlue<=0) = 0;


handles.aoLaserSession.outputSingleScan(voltPowerBlue);
guidata(hObject, handles);



function red_laser_output_power_Callback(hObject, eventdata, handles)
inPreview = get(handles.preview, 'Value');

%In preview mode, laser power is set to 2 if not set to 1 or 0
if (inPreview == 0)
    laserPowerBlue = str2num(get(handles.blue_laser_output_power, 'String'));
    laserPowerRed = str2num(get(handles.red_laser_output_power, 'String'));
else
    if(str2num(get(handles.blue_laser_output_power, 'String'))>=2)
        laserPowerBlue = 2;
    else
        laserPowerBlue = 0;
    end
    if(str2num(get(handles.red_laser_output_power, 'String'))>=2)
        laserPowerRed = 2;
    else
        laserPowerRed = 0;
    end
end

% Old formula (before changing to 2000Ohm input impedance)
% voltPower = laserPower*0.0694-0.1140;
% NEW FORMULA:)
voltPowerBlue =(laserPowerBlue*laserPowerBlue)*0.0001+(laserPowerBlue*0.0305)+0.0012;  %SCAPE4 lasers
%voltPowerBlue =(laserPowerBlue*laserPowerBlue)*1.423e-05+(laserPowerBlue*0.03241)+0.002877; %Uptown SCAPE lasers
voltPowerBlue(voltPowerBlue<=0) = 0;

%voltPowerRed =
%(laserPowerRed*laserPowerRed)*0.0001+(laserPowerRed*0.0349)+0.0107; %SCAPE4 lasers
voltPowerRed = (laserPowerRed*laserPowerRed)*3.08e-05+(laserPowerRed*0.03788)-0.003457; %Uptown SCAPE lasers
voltPowerRed(voltPowerRed<=0) = 0;
%if(voltPowerRed>2)
%     voltPowerRed = 2;
%     set(handles.red_laser_output_power, 'String', '47.12');
%end

%KPedit - if both lasers are checked,and only one AO is used for both, the voltage assigned to the blue laser will be used for both
if (get(handles.blue_laser_checkbox, 'Value')==0)
    handles.aoLaserSession.outputSingleScan(voltPowerRed);
else
    handles.aoLaserSession.outputSingleScan(voltPowerBlue);
end

guidata(hObject, handles);

function step_size_Callback(hObject, eventdata, handles)
handles.stepSize = str2num(get(handles.step_size, 'String'));
scanFOV = handles.scanFOV;
framerate = str2num(get(handles.camera_framerate, 'String'));

% Makes sure the FOV is evenly divisible by the step size and updates the
% FOV if this isn't the case
temp = ceil(scanFOV/handles.stepSize);
handles.scanFOV = handles.stepSize*temp;
set(handles.scan_FOV, 'String', num2str(handles.scanFOV));

% Updates the number of scan steps this will take
handles.numScanSteps = temp+1;
set(handles.num_scan_steps, 'String', num2str(handles.numScanSteps));

% Adds 1 for flyback
numFramesPerVolume = handles.numScanSteps+1;

% Updates the scan rate
handles.scanRate = framerate/numFramesPerVolume;
set(handles.volumetric_scan_rate, 'String', num2str(handles.scanRate));
handles.paramset = 1;
scan_calculations(hObject, handles);



function step_size2_Callback(hObject, eventdata, handles)
% hObject    handle to step_size2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.stepSize2 = str2num(get(handles.step_size2, 'String'));
scanFOV2 = handles.scanFOV2;
framerate2 = str2num(get(handles.camera_framerate2, 'String'));

% Makes sure the FOV is evenly divisible by the step size and updates the
% FOV if this isn't the case
temp2 = ceil(scanFOV2/handles.stepSize2);
handles.scanFOV2 = handles.stepSize2*temp2;
set(handles.scan_FOV2, 'String', num2str(handles.scanFOV2));

% Updates the number of scan steps this will take
handles.numScanSteps2 = temp2+1;
set(handles.num_scan_steps2, 'String', num2str(handles.numScanSteps2));

% Adds 1 for flyback
numFramesPerVolume2 = handles.numScanSteps2+1;

% Updates the scan rate
handles.scanRate2 = framerate2/numFramesPerVolume2;
set(handles.volumetric_scan_rate2, 'String', num2str(handles.scanRate2));
handles.paramset = 1;
scan_calculations(hObject, handles);



function volumetric_scan_rate_Callback(hObject, eventdata, handles)
handles.scanRate = str2num(get(handles.volumetric_scan_rate, 'String'));
framerate   = str2num(get(handles.camera_framerate, 'String'));
scanFOV     = handles.scanFOV;
numFramesPerVolume = floor(framerate/handles.scanRate);

% Find/update the closest possible scan rate
newScanRate = framerate/numFramesPerVolume;
handles.scanRate = newScanRate;
set(handles.volumetric_scan_rate, 'String', num2str(handles.scanRate));

% Update num scan steps
handles.numScanSteps = numFramesPerVolume-1;
set(handles.num_scan_steps, 'String', num2str(handles.numScanSteps));

% Update step size
handles.stepSize = scanFOV/(handles.numScanSteps-1);
set(handles.step_size, 'String', num2str(handles.stepSize));
handles.paramset = 1;
scan_calculations(hObject, handles);



function volumetric_scan_rate2_Callback(hObject, eventdata, handles)
% hObject    handle to volumetric_scan_rate2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.scanRate2 = str2num(get(handles.volumetric_scan_rate2, 'String'));
framerate2   = str2num(get(handles.camera_framerate2, 'String'));
scanFOV2     = handles.scanFOV2;
numFramesPerVolume2 = floor(framerate2/handles.scanRate2);

% Find/update the closest possible scan rate
newScanRate2 = framerate2/numFramesPerVolume2;
handles.scanRate2 = newScanRate2;
set(handles.volumetric_scan_rate2, 'String', num2str(handles.scanRate2));

% Update num scan steps
handles.numScanSteps2 = numFramesPerVolume2-1;
set(handles.num_scan_steps2, 'String', num2str(handles.numScanSteps2));

% Update step size
handles.stepSize2 = scanFOV2/(handles.numScanSteps2-1);
set(handles.step_size2, 'String', num2str(handles.stepSize2));
handles.paramset = 1;
scan_calculations(hObject, handles);


function num_scan_steps_Callback(hObject, eventdata, handles)
handles.numScanSteps        = str2double(get(handles.num_scan_steps, 'String'));
cameraFrameRate     = str2num(get(handles.camera_framerate, 'String'));
scanFOV             = handles.scanFOV;

% Update step size
handles.stepSize = scanFOV/(handles.numScanSteps-1);
set(handles.step_size, 'String', num2str(handles.stepSize));

% Extra step for the flyback scan
numFramesPerVolume = handles.numScanSteps+1;
% Calculate new volume rate and update
newVolumeRate = cameraFrameRate/numFramesPerVolume;
handles.scanRate = newVolumeRate;
set(handles.volumetric_scan_rate, 'String', num2str(newVolumeRate));
handles.paramset = 1;
scan_calculations(hObject, handles);



function num_scan_steps2_Callback(hObject, eventdata, handles)
% hObject    handle to num_scan_steps2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.numScanSteps2        = str2double(get(handles.num_scan_steps2, 'String'));
cameraFrameRate2     = str2num(get(handles.camera_framerate2, 'String'));
scanFOV2             = handles.scanFOV2;

% Update step size
handles.stepSize2 = scanFOV2/(handles.numScanSteps2-1);
set(handles.step_size2, 'String', num2str(handles.stepSize2));

% Extra step for the flyback scan
numFramesPerVolume2 = handles.numScanSteps2+1;
% Calculate new volume rate and update
newVolumeRate2 = cameraFrameRate2/numFramesPerVolume2;
handles.scanRate2 = newVolumeRate2;
set(handles.volumetric_scan_rate2, 'String', num2str(newVolumeRate2));
handles.paramset = 1;
scan_calculations(hObject, handles);




function camera_framerate_Callback(hObject, eventdata, handles)

% Update camera frame rate
handles.cameraFramerate = str2num(get(handles.camera_framerate, 'String'));
handles.framerate_newInput = 1;
numScanSteps            = handles.numScanSteps;

% Add frame for flyback
numFramesPerVolume = numScanSteps+1;
% Calculate new volume rate and store
newVolumeRate = handles.cameraFramerate/numFramesPerVolume;
% Update parameters
handles.scanRate = newVolumeRate;
set(handles.volumetric_scan_rate, 'String', num2str(newVolumeRate));
%if get(handles.max_frame_rate,'Value')
%Marking framerate to be recalculated
handles.framerate_newFPS = 1;
%end
handles.paramset = 1;
scan_calculations(hObject, handles);

function camera_framerate2_Callback(hObject, eventdata, handles)

% Update camera frame rate
handles.cameraFramerate2 = str2num(get(handles.camera_framerate2, 'String'));
handles.framerate_newInput = 1;
numScanSteps2            = handles.numScanSteps2;

% Add frame for flyback
numFramesPerVolume2 = numScanSteps2+1;
% Calculate new volume rate and store
newVolumeRate2 = handles.cameraFramerate2/numFramesPerVolume2;
% Update parameters
handles.scanRate2 = newVolumeRate2;
set(handles.volumetric_scan_rate2, 'String', num2str(newVolumeRate2));
%Untick max frame rate box
set(handles.max_frame_rate, 'Value', 0);
%if get(handles.max_frame_rate,'Value')
%Marking framerate to be recalculated
handles.framerate_newFPS = 1;
%end
handles.paramset = 1;
scan_calculations(hObject, handles);


function scan_length_Callback(hObject, eventdata, handles)
handles.scanLength = str2num(get(handles.scan_length, 'String'));
handles.paramset = 1;
scan_calculations(hObject, handles);



function scan_length2_Callback(hObject, eventdata, handles)
handles.scanLength2 = str2num(get(handles.scan_length2, 'String'));
handles.paramset = 1;
scan_calculations(hObject, handles);



function frame_rate_changed(hObject, eventdata, handles)
% Update camera frame rate
fileName = [handles.GUI_filepath, 'cameraConfig_frameRateUpdate.txt'];
fid = fopen(fileName);

isHR = strcmp('HR', handles.tgroup.SelectedTab.Title);
if(isHR == 1)
    handles.cameraFramerate = str2num(strtrim(fgets(fid)));
    set(handles.camera_framerate, 'String', num2str(handles.cameraFramerate));
    fclose(fid);
    numScanSteps            = str2num(get(handles.num_scan_steps, 'String'));
    % Add frame for flyback
    numFramesPerVolume = numScanSteps+1;
    % Update parameters
    handles.scanRate = handles.cameraFramerate/numFramesPerVolume;
    set(handles.volumetric_scan_rate, 'String', num2str(handles.scanRate));
    
else
    handles.cameraFramerate2 = str2num(strtrim(fgets(fid)));
    set(handles.camera_framerate2, 'String', num2str(handles.cameraFramerate2));
    fclose(fid);
    numScanSteps            = str2num(get(handles.num_scan_steps2, 'String'));
    % Add frame for flyback
    numFramesPerVolume = numScanSteps+1;
    % Update parameters
    handles.scanRate2 = handles.cameraFramerate2/numFramesPerVolume;
    set(handles.volumetric_scan_rate2, 'String', num2str(handles.scanRate2));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Galvo Related Functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scan_FOV_Callback(hObject, eventdata, handles)
handles.scanFOV = str2double(get(handles.scan_FOV, 'String'));
offset = str2double(get(handles.galvo_um_offset, 'String'));
minPosition = offset-(handles.scanFOV/2);
maxPosition = offset+(handles.scanFOV/2);

% Update step size
numScanSteps = str2double(get(handles.num_scan_steps, 'String'));

% I remind you again, dear coder, [0 1 2 3 4 5] spans a range of 5, has a step size of 1, and
% includes 6 elements
handles.stepSize = handles.scanFOV/(numScanSteps-1);
set(handles.step_size, 'String', num2str(handles.stepSize));

% Update slider labels
set(handles.galvo_position_max, 'String', [num2str(maxPosition) ' um']);
set(handles.galvo_position_min, 'String', [num2str(minPosition) ' um']);

% Update slider limits and set position to offset
set(handles.galvo_position_slider, 'Max', maxPosition);
set(handles.galvo_position_slider, 'Min', minPosition);

guidata(hObject, handles);

% Will check whether the current position is out of bounds and if so,
% adjust the current position accordingly. Otherwise, nothing will happen
shutterValue = get(handles.blue_shutter, 'Value');
if (shutterValue == 0)
    galvo_position_display_Callback(hObject, eventdata, handles);
else
    galvo_position_display_Callback(hObject, eventdata, handles);
end
handles.paramset = 1;
scan_calculations(hObject, handles);




function scan_FOV2_Callback(hObject, eventdata, handles)
% hObject    handle to scan_FOV2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scan_FOV as text
%        str2double(get(hObject,'String')) returns contents of scan_FOV as a double
handles.scanFOV2 = str2double(get(handles.scan_FOV2, 'String'));
offset = str2double(get(handles.galvo_um_offset, 'String'));
minPosition2 = offset-(handles.scanFOV2/2);
maxPosition2 = offset+(handles.scanFOV2/2);

% Update step size
numScanSteps2 = str2double(get(handles.num_scan_steps2, 'String'));

% I remind you again, dear coder, [0 1 2 3 4 5] spans a range of 5, has a step size of 1, and
% includes 6 elements
handles.stepSize2 = handles.scanFOV2/(numScanSteps2-1);
set(handles.step_size2, 'String', num2str(handles.stepSize2));

% Update slider labels
set(handles.galvo_position_max, 'String', [num2str(maxPosition2) ' um']);
set(handles.galvo_position_min, 'String', [num2str(minPosition2) ' um']);

% Update slider limits and set position to offset
set(handles.galvo_position_slider, 'Max', maxPosition2);
set(handles.galvo_position_slider, 'Min', minPosition2);

guidata(hObject, handles);

% Will check whether the current position is out of bounds and if so,
% adjust the current position accordingly. Otherwise, nothing will happen
shutterValue = get(handles.blue_shutter, 'Value');
if (shutterValue == 0)
    galvo_position_display_Callback(hObject, eventdata, handles);
else
    galvo_position_display_Callback(hObject, eventdata, handles);
end
handles.paramset = 1;
scan_calculations(hObject, handles);


function galvo_position_display_Callback(hObject, eventdata, handles)

% Get new position
newpos = str2double( get(handles.galvo_position_display, 'String'));
maxPosition = get(handles.galvo_position_slider, 'Max');
minPosition = get(handles.galvo_position_slider, 'Min');
offset = str2double(get(handles.galvo_um_offset, 'String'));
if ((newpos > maxPosition) || (newpos < minPosition))
    % Set new position to max or min value
    if (newpos > offset)
        newpos = maxPosition;
    else
        newpos = minPosition;
    end
    
    % Update position display
    set(handles.galvo_position_display, 'String', num2str(newpos, 3));
end

% Update slider with new position value
set(handles.galvo_position_slider, 'Value', newpos);

% Run slider callback to reposition galvo
galvo_position_slider_Callback(hObject, eventdata, handles);

function center_galvo_Callback(hObject, eventdata, handles)

offset = get(handles.galvo_um_offset, 'String');
set(handles.galvo_position_display, 'String', offset);

galvo_position_display_Callback(hObject, eventdata, handles);

function galvo_um_offset_Callback(hObject, eventdata, handles)
isHR = strcmp('HR', handles.tgroup.SelectedTab.Title);
if(isHR == 1)
    handles.scanFOV = str2double(get(handles.scan_FOV, 'String'));
    offset = str2double(get(handles.galvo_um_offset, 'String'));
    minPosition = offset-(handles.scanFOV/2);
    maxPosition = offset+(handles.scanFOV/2);
else
    handles.scanFOV2 = str2double(get(handles.scan_FOV2, 'String'));
    offset = str2double(get(handles.galvo_um_offset, 'String'));
    minPosition = offset-(handles.scanFOV2/2);
    maxPosition = offset+(handles.scanFOV2/2);
end


% Update slider labels
set(handles.galvo_position_max, 'String', [num2str(maxPosition) ' um']);
set(handles.galvo_position_min, 'String', [num2str(minPosition) ' um']);

% Update slider limits and set position to offset
set(handles.galvo_position_slider, 'Max', maxPosition);
set(handles.galvo_position_slider, 'Min', minPosition);

% update position display properties
set(handles.galvo_position_display, 'String', offset);

guidata(hObject, handles);

% Call galvo position display. This will update that box and the slider
% value
galvo_position_display_Callback(hObject, eventdata, handles);

scan_calculations(hObject,handles);


function galvo_position_slider_Callback(hObject, eventdata, handles)
%h = waitbar(0, 'Moving to new position please wait.');  % okay - don't
%reallyyyy need this..

% get new position
conversionFactor = handles.info.GUIcalFactors.xK_umPerVolt;
newpos = get(handles.galvo_position_slider, 'Value');
globalOffset = handles.globalGalvoOffset; % Value in microns (set by me when I install the system)

% move to position
% outputVoltage = [0 (newpos+globalOffset)/conversionFactor voltPower];
% kpedit - no stim on uptown system
%outputVoltage = [0 (newpos+globalOffset)/conversionFactor];
outputVoltage = (newpos+globalOffset)/conversionFactor;
handles.aoSession.outputSingleScan(outputVoltage);

% update position display
handles.currentGalvoPosition = newpos;
guidata(hObject, handles);
set(handles.galvo_position_display, 'String', num2str(handles.currentGalvoPosition, 3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Pure Code Functions %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function stop_scan_Callback(hObject, eventdata, handles)
% close shutters
handles.dioSession.outputSingleScan([0 0 0]);

% stop AI
if (1 == handles.aiSession.IsRunning)
    handles.aiSession.stop;
    disp('AI stopped!');
end

% stop AO
if (1 == handles.aoSession.IsRunning)
    handles.aoSession.stop;
    disp('AO stopped!');
end

%reset solis
system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_stopPreview.ahk');
pause(1)
system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_setParameters.ahk');
pause(2)
system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\returnToGUI.ahk');

% stop and restore normal function to motorized stage and update positions
isHR = strcmp('HR', handles.tgroup.SelectedTab.Title);
if (get(handles.motorstage, 'Value') && isHR==0)
    handles.stage.SetKCubeTriggerParams(0, 0, 1, 0, 1); %will move to abs position without trigger
    handles.stage.EnableHWChannel(0);
    handles.stage.StopImmediate(0); %stop motion if not already stopped
    % Verify stage position
    pause(1)
    pos = handles.stage.GetPosition_Position(0);
    set(handles.stage_slider, 'Value', pos);
    set(handles.stage_pos, 'String', num2str(pos));  
end

function ai_channels_to_plot_Callback(hObject, eventdata, handles)
if (1 == isfield(handles, 'lastStimData'))
    channelsToPlot = 1 + str2num( get(handles.ai_channels_to_plot, 'String'));
    plot(handles.stim_axes, handles.lastStimData.time, handles.lastStimData.channels(channelsToPlot, :));
    
    legendText = [];
    if (find(channelsToPlot == 1))
        legendText = {'Galvo Feedback (ignore)'};
    end
    if (find(channelsToPlot == 2))
        legendText = [legendText {'Galvo Drive'}];
    end
    if (find(channelsToPlot == 3))
        legendText = [legendText {'Camera ds'}];
    end
    
    legend(handles.stim_axes, legendText, 'Location', 'NorthEast');
end

function writeConfigFile(hObject, handles)
presentDirectory = pwd;
cd(handles.GUI_filepath);
fileName = 'cameraConfig.txt';
% Want to read certain lines and change the value of those lines

fid = fopen(fileName, 'w');

% Write in the following order:
% fps
% width
% height
% left
% binning value
% scan duration (sec)
% directory
% scanName

info = handles.info;
dataDirectory = info.dataDirectory;
indices = find(dataDirectory == '\');
isHR = strcmp('HR', handles.tgroup.SelectedTab.Title);
isFullFrame = get(handles.fullFrame, 'Value');
inPreview = get(handles.preview, 'Value');

for i = 1:length(indices)-1
    dataDirectory = [dataDirectory(1:indices(i)) '\' dataDirectory(indices(i)+1:end)];
    indices = indices+1;
end
dataDirectory = [dataDirectory '\'];

if (handles.info.loadParameters.maxFrameRate == 1 && isHR == 0 && isFullFrame == 0 && inPreview == 0)
    framerate = 3000;
else
    framerate = handles.info.loadParameters.desiredCameraFramerate;
end
% % % saves the current scan's time stamp to the info file
info.scanStartTimeApprox = datestr(now);

fprintf(fid, [num2str(framerate) '\n']);
fprintf(fid, [num2str(info.camera.xROI) '\n']);
fprintf(fid, [num2str(info.camera.yROI) '\n']);
fprintf(fid, [num2str(info.camera.x_left) '\n']);
fprintf(fid, [num2str(info.camera.binFactor) '\n']);
fprintf(fid, [num2str(info.daq.scanLength) '\n']);
fprintf(fid, [dataDirectory '\n']);
fprintf(fid, [info.scanName '\n']);

fclose(fid);

cd(presentDirectory);


function blue_shutter_Callback(hObject, eventdata, handles)
blueShutterValue = get(handles.blue_shutter, 'Value');
redShutterValue = get(handles.red_shutter, 'Value');

if (blueShutterValue == 1)
    disp('Blue Shutter Open');
else
    disp('Blue Shutter Closed');
end
digitalLines = [0 blueShutterValue redShutterValue]; % Set blue channel to high if button pressed
handles.dioSession.outputSingleScan(digitalLines);

function red_shutter_Callback(hObject, eventdata, handles)
blueShutterValue = get(handles.blue_shutter, 'Value');
redShutterValue = get(handles.red_shutter, 'Value');

if (redShutterValue == 1)
    disp('Red Shutter Open');
else
    disp('Red Shutter Closed');
end

digitalLines = [0 blueShutterValue redShutterValue]; % Set blue channel to high if button pressed
handles.dioSession.outputSingleScan(digitalLines);

function update_extra_info_parameters(hObject, handles)
% This function saves parameters to an info file that 1) aren't used for
% data loading and 2) are required to load the info file for future
% acquisitions
isHR = strcmp('HR', handles.tgroup.SelectedTab.Title);
if(isHR == 1)
    handles.cameraFramerate = str2num(get(handles.camera_framerate, 'String'));
    handles.stepSize = str2num(get(handles.step_size, 'String'));
    handles.scanLength = str2num(get(handles.scan_length, 'String'));
else
    handles.cameraFramerate2 = str2num(get(handles.camera_framerate2, 'String'));
    handles.stepSize2 = str2num(get(handles.step_size2, 'String'));
    handles.scanLength2 = str2num(get(handles.scan_length2, 'String'));
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Stim Related %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Scan Calculations callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function experiment_notes_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function objective_list_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);


function pixel_binning_Callback(hObject, eventdata, handles)
handles.paramset = 1;
scan_calculations(hObject, handles);

function pixel_binning2_Callback(hObject, eventdata, handles)
handles.paramset = 1;
scan_calculations(hObject, handles);

function scan_name_stem_Callback(hObject, eventdata, handles)
handles.paramset = 1;
scan_calculations(hObject, handles);

function data_directory_Callback(hObject, eventdata, handles)
handles.paramset = 1;
scan_calculations(hObject, handles);

function x_ROI_Callback(hObject, eventdata, handles)
set(handles.x_ROI2, 'String', get(handles.x_ROI, 'String'))
handles.paramset = 1;
scan_calculations(hObject, handles);

function y_ROI_Callback(hObject, eventdata, handles)
handles.paramset = 1;
scan_calculations(hObject, handles);

function x_ROI2_Callback(hObject, eventdata, handles)
set(handles.x_ROI, 'String', get(handles.x_ROI2, 'String'))
handles.paramset = 1;
scan_calculations(hObject, handles);

function y_ROI2_Callback(hObject, eventdata, handles)
if get(handles.max_frame_rate,'Value')
%Marking framerate to be recalculated
handles.framerate_newFPS = 1;
end
handles.paramset = 1;
scan_calculations(hObject, handles);

function camera_binning_Callback(hObject, eventdata, handles)
handles.paramset = 1;
scan_calculations(hObject, handles);


function ai_input_channels_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function ai_sample_rate_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function HR_checkbox_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

% --- Executes on button press in max_frame_rate.
function max_frame_rate_Callback(hObject, eventdata, handles)
if get(handles.max_frame_rate,'Value')
%Marking framerate to be recalculated
handles.framerate_newFPS = 1;
end
handles.paramset = 1;
scan_calculations(hObject, handles)

function stimTable_CellEditCallback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function camera_left_Callback(hObject, eventdata, handles)
set(handles.camera_left2, 'String', get(handles.camera_left, 'String'))
handles.paramset = 1;
scan_calculations(hObject, handles);

function camera_left2_Callback(hObject, eventdata, handles)
set(handles.camera_left, 'String', get(handles.camera_left2, 'String'))
handles.paramset = 1;
scan_calculations(hObject, handles);

function optical_scan_angle_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function blue_laser_checkbox_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function red_laser_checkbox_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kineticSeriesLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kineticSeriesLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function camera_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camera_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lateral_FOV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lateral_FOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function depth_FOV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to depth_FOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function scan_FOV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scan_FOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function num_scan_steps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_scan_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function volumes_to_collect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to volumes_to_collect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter_wheel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_wheel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function experiment_notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to experiment_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pixel_binning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixel_binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function camera_framerate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camera_framerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ai_channels_to_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ai_channels_to_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOT NEEDED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pre_stim_time_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function stim_time_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function post_stim_time_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function stim_rate_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function stim_width_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function use_stim_checkbox_Callback(hObject, eventdata, handles)
useStim = get(handles.use_stim_checkbox, 'Value');
if(useStim == 1)
    set(handles.stimPanel, 'Visible', 'on');
    set(handles.analogdatainput, 'Visible', 'on');
    %     set(findall(handles.stimPanel, '-property', 'enable'), 'enable', 'on');
else
    set(handles.stimPanel, 'Visible', 'off');
    set(handles.analogdatainput, 'Visible', 'off');
    %     set(findall(handles.stimPanel, '-property', 'enable'), 'enable', 'off');
end
scan_calculations(hObject, handles);

function stim_sample_rate_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function two_channel_checkbox_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function sawtooth_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function display_scanner_waveform_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function camera_y_bottom_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function alignment_checkBox_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function depth_FOV_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function lateral_FOV_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function kineticSeriesLength_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function set_fastest_rowVal_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function green_laser_checkbox_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);


function pixels_per_line_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);


function edit51_Callback(hObject, eventdata, handles)% --> what is this?
% hObject    handle to camera_framerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of camera_framerate as text
%        str2double(get(hObject,'String')) returns contents of camera_framerate as a double

function volumes_to_collect_Callback(hObject, eventdata, handles)
isHR = strcmp('HR', handles.tgroup.SelectedTab.Title);
if(isHR == 1)
    volumetricScanRate = str2double(get(handles.volumetric_scan_rate, 'String'));
    numVolumesAcquired = str2double(get(handles.volumes_to_collect, 'String'));
    scanLength = numVolumesAcquired/volumetricScanRate;
    set(handles.scan_length, 'String', num2str(scanLength));
else
    volumetricScanRate = str2double(get(handles.volumetric_scan_rate2, 'String'));
    numVolumesAcquired = str2double(get(handles.volumes_to_collect, 'String'));
    scanLength = numVolumesAcquired/volumetricScanRate;
    set(handles.scan_length2, 'String', num2str(scanLength));
end
scan_calculations(hObject, handles);


function green_shutter_Callback(hObject, eventdata, handles)
shutterValue = get(handles.green_shutter, 'Value');
if (1 == shutterValue)
    handles.info.shutterToggle(2) = 1;
    disp('Green Shutter Open');
elseif (0 == shutterValue)
    handles.info.shutterToggle(2) = 0;
    disp('Green Shutter Open');
end
guidata(hObject, handles);

handles.dioSession.outputSingleScan([0 handles.info.shutterToggle]);



function edit52_Callback(hObject, eventdata, handles)
% hObject    handle to ai_channels_to_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ai_channels_to_plot as text
%        str2double(get(hObject,'String')) returns contents of ai_channels_to_plot as a double

function loadAlignment(hObject, eventdata, handles)
dataDirectory = handles.info.dataDirectory;
scanName = handles.info.scanName;
infoFilePath = fullfile(dataDirectory, [scanName '_info.mat']);

% cd (dataDirectory)
if exist(infoFilePath)
    % load info file
    disp('LOADING INFO FILE')
    load(infoFilePath);
    disp(info.daq)
    
    % Scake data set size by spatial bin factor
    binFactor = info.camera.binFactor;          % Camera Bin Number
    numDepths = info.camera.yROI / binFactor;   % Number of Depth Pixels
    numLatPix = info.camera.xROI / binFactor;     % Number of Lateral Pixels
    %         numDepths = info.camera.yROI;
    %         numLatPix = info.camera.xROI;
    scanRate = info.daq.scanRate;               % Volumetric Scan Rate (VPS)
    
    
    % Given in microns/pixel for depth, lateral and scan dimensions (in that
    % order)
    % Set the conversion factors based on objective used
    presentDirectory = pwd;
    %     cd C:\Users\Hodor\Documents\LSIPT_galvo_drive_GUI;
    load savedConversionFactors.mat
    cd(presentDirectory)
    objective = info.objective;
    
    
    if (isempty(strfind(objective, '10x')) == 0)
        cFact = savedConversionFactors(1);
        conversionFactors(1) = cFact.depth_umPerPix;     % Does change as a function of objective magnification
        conversionFactors(2 ) = cFact.lateral_umPerPix;   % Does change as a function of objective magnification
        conversionFactors(3) = 1/(cFact.scan_VoltsPerUm*info.daq.pixelsPerLine/info.daq.scanAngle); % Doesn't change as a function of objective magnification (Volts to galvo/um)alvo/um)
    elseif(isempty(strfind(objective, '20x')) == 0)
        cFact = savedConversionFactors(2);
        conversionFactors(1) = cFact.depth_umPerPix;     % Does change as a function of objective magnification
        conversionFactors(2) = cFact.lateral_umPerPix;   % Does change as a function of objective magnification
        conversionFactors(3) = 1/(cFact.scan_VoltsPerUm*info.daq.pixelsPerLine/info.daq.scanAngle); % Doesn't change as a function of objective magnification (Volts to galvo/um)
    end
    
    
    %conversionFactors = 1./[0.4724/binFactor .5512/binFactor 0.0028*info.daq.pixelsPerLine/info.daq.scanAngle];
    %conversionFactors = 1./[1.414/binFactor 1.0709/binFactor 0.0031*info.daq.pixelsPerLine/info.daq.scanAngle];
    
    disp('LOADING METADATA')
    % Creates file path of zyla metadata file (.ini)
    zylaInfoFilePath = fullfile(dataDirectory, scanName, 'acquisitionmetadata.ini');
    FID = fopen(zylaInfoFilePath, 'r');
    zylaMetaData = fread(FID, '*char')';
    fclose(FID);
    % Reads the number of frames in a spool file
    startIndex = strfind(zylaMetaData, 'ImagesPerFile') + length('ImagesPerFile = ');
    numFramesPerSpool = str2double(zylaMetaData(startIndex:end));
    % Reads the image size in the spool file
    startIndex = strfind(zylaMetaData, 'ImageSize') + length('ImageSize = ');
    ImageSize = str2double(zylaMetaData(startIndex : startIndex + 7));
    
    numFrames = info.daq.pixelsPerLine;  % Total number of frames acquired
    numFramesAcquired = round(str2num(info.camera.kineticSeriesLength));
    
    % add extra columns for 2 buffer rows (why?)
    if (2 == info.camera.binFactor)
        numColumns = numDepths + 1;
        numRows = ImageSize / 2 / numColumns;
        if(mod(numRows, 1) ~= 0)
            numColumns = numDepths+2;
            numRows = ImageSize/2/numColumns;
        end
    else
        numColumns = numDepths + 2;
        numRows = numLatPix;
    end
    % Creates a list of spool file names that exist for a given run
    disp('CREATING SPOOL FILE NAMES')
    spoolCounter=0;
    for i = 1:floor(numFramesAcquired/(numFramesPerSpool))
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
    filesToLoad = namesOut
    clear namesOut tempName temp a spoolCounter
    
    cd(dataDirectory);
    cd(scanName);
    disp('ALLOCATING MEMORY')
    tic
    moviesub = uint16(0);
    moviesub(numRows, numColumns, numFramesAcquired) = 0;
    toc
    spoolFileSize = numRows*numColumns*numFramesPerSpool;
    exitCondition = 0;
    fileCounter = 1;
    h = waitbar(0, 'Please Wait');
    movieSize = prod(size(moviesub));
    
    while exitCondition == 0
        
        firstIndex = (fileCounter-1)*spoolFileSize +1;
        lastIndex = (fileCounter*spoolFileSize);
        if (firstIndex > movieSize) % All of these pixels will be dark images
            exitCondition = 1;
        elseif (lastIndex > movieSize) % Some of these pixels will be dark images
            dataLeftToFill = [firstIndex:movieSize];
            if (length(dataLeftToFill) ~= 0)
                fileToRead = filesToLoad{fileCounter};
                fid = fopen(fileToRead, 'r');
                X = fread(fid, 'uint16=>uint16')';
                fclose(fid);
                moviesub(dataLeftToFill) = ...
                    X(1:length(dataLeftToFill));
            end
            exitCondition = 1;
        else % All of these pixels will be moviesub
            fileToRead = filesToLoad{fileCounter};
            fid = fopen(fileToRead, 'r');
            moviesub(firstIndex:lastIndex) = ...
                fread(fid, 'uint16=>uint16')';
            fclose(fid);
            
        end
        waitbar(lastIndex/movieSize)
        fileCounter = fileCounter +1;
    end
    
    
    darkImage = moviesub(:, :, numFrames+1:end);
    meanDark = squeeze(mean(mean(darkImage, 1), 2));
    toKeep = find(meanDark < 110);
    darkImage = darkImage(:, :, toKeep);
    moviesub(:, :, numFrames+1:end) = [];
    
    close (h)
else
    keyboard
end

FOV = size(moviesub).*[conversionFactors(2) conversionFactors(1) conversionFactors(3)];
moviesub(:, end-1:end, :, :) = [];
darkImage(:, end-1:end, :, :) = [];
cd(dataDirectory)

figure(100);
subplot(131)
imagesc(squeeze(max(moviesub, [], 1)));
subplot(132)
imagesc(squeeze(max(moviesub, [], 2)));
subplot(133)
imagesc(squeeze(max(moviesub, [], 3)));

clear moviesub darkImage
disp('FINISHED LOADING IN DATA');
%close all

% --- Executes during object creation, after setting all properties.
function galvo_um_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to galvo_um_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function laser_output_power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_laser_output_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function red_laser_output_power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_laser_output_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function blue_laser_output_power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_laser_output_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in preview.
function preview_Callback(hObject, eventdata, handles)
% hObject    handle to preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
%set laser power to 2 unless power is already 0 (both lasers)
%sets framerate to 10
%restores previous laser power/framerate upon unclicking
% 2 main modifications: laser power and frame rate
isFullFrame = get(handles.fullFrame, 'Value');
inPreview = get(handles.preview, 'Value');
isHR = strcmp('HR', handles.tgroup.SelectedTab.Title);

if(inPreview == 1)
    if (get(handles.solispreview, 'Value') == 0)
        set(handles.solispreview, 'Value', 1)
    else
        system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_stopPreview.ahk');
        system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\returnToGUI.ahk');
    end
   
    %Disable fullframe mode
    if(isFullFrame == 1)
        set(handles.fullFrame, 'Value', 0);
        if(isHR == 1)
            set(handles.x_ROI, 'Enable', 'on');
            set(handles.x_ROI, 'BackgroundColor', [1 1 1]);
            set(handles.y_ROI, 'Enable', 'on');
            set(handles.y_ROI, 'BackgroundColor', [1 1 1]);
            set(handles.camera_left, 'Enable', 'on');
            set(handles.camera_left, 'BackgroundColor', [1 1 1]);
            %Enabling framerate field, max framerate box
            set(handles.camera_framerate, 'Enable', 'on');
            set(handles.camera_framerate, 'BackgroundColor', [1 1 1]);
        else
            set(handles.x_ROI2, 'Enable', 'on');
            set(handles.x_ROI2, 'BackgroundColor', [1 1 1]);
            set(handles.y_ROI2, 'Enable', 'on');
            set(handles.y_ROI2, 'BackgroundColor', [1 1 1]);
            set(handles.camera_left2, 'Enable', 'on');
            set(handles.camera_left2, 'BackgroundColor', [1 1 1]);
            %Enabling framerate field, max framerate box
            set(handles.camera_framerate2, 'Enable', 'on');
            set(handles.camera_framerate2, 'BackgroundColor', [1 1 1]);
            set(handles.max_frame_rate, 'Enable', 'on');
            set(handles.max_frame_rate, 'BackgroundColor', [1 1 1]);
        end
        previewwason = 1;
    else
        previewwason = 0;
    end
    
    %Disable laser power fields
    set(handles.blue_laser_output_power, 'Enable', 'off');
    set(handles.red_laser_output_power, 'Enable', 'off');
    %Disable framerate field and max framerate box
    if(isHR==1)
        set(handles.camera_framerate, 'Enable', 'off');
        if ~previewwason
            handles.originalFPS = handles.cameraFramerate;
        end
    else
        set(handles.camera_framerate2, 'Enable', 'off');
        set(handles.max_frame_rate, 'Enable', 'off');
        if ~previewwason
        handles.originalFPS2 = handles.cameraFramerate2;
        end
    end
    %set filter wheel to 1%
%     handles.originalfiltval = get(handles.filter_wheel, 'Value');
%      if get(handles.highpower, 'Value')
%         set(handles.filter_wheel, 'Value', 3);
%         input = 'pos=3';
%     elseif get(handles.medpower, 'Value')
%         set(handles.filter_wheel, 'Value', 3);
%         input = 'pos=3';
%     else
%         set(handles.filter_wheel, 'Value', 4);
%         input = 'pos=4';
%     end
%     fprintf(handles.filterWheelObj, input);

userChoice = num2str(get(handles.alignpower, 'Value'));
input = ['pos=' userChoice];
fprintf(handles.filterWheelObj, input);
    
else
    if (get(handles.solispreview, 'Value') == 1)
        set(handles.solispreview, 'Value', 0)
        system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_stopPreview.ahk');
        system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\returnToGUI.ahk');
    end  
    
    %Enable laser power fields
    set(handles.blue_laser_output_power, 'Enable', 'on');
    set(handles.red_laser_output_power, 'Enable', 'on');
    %Enable framerate field and max framerate box.
    if (isHR==1)
        set(handles.camera_framerate, 'Enable', 'on');
        handles.cameraFramerate = handles.originalFPS;
        set(handles.camera_framerate, 'String', num2str(handles.originalFPS));
    else
        set(handles.camera_framerate2, 'Enable', 'on');
        set(handles.max_frame_rate, 'Enable', 'on');
        handles.cameraFramerate2 = handles.originalFPS2;
        set(handles.camera_framerate2, 'String', num2str(handles.originalFPS2));
    end
    %Marking framerate to be recalculated
    handles.framerate_newInput = 1;
    %reset filter wheel 
    set(handles.filter_wheel, 'Value', handles.originalfiltval);
    input = ['pos=' num2str(handles.originalfiltval)];
    fprintf(handles.filterWheelObj, input);
end

% Update handles structure
guidata(hObject, handles);

%reset laser power
blue_laser_output_power_Callback(hObject, eventdata, handles);
red_laser_output_power_Callback(hObject, eventdata, handles);

%turn on red and blue lasers if alignment mode is turned on, (and off if turned off) 
if inPreview == 1
    if get(handles.blue_laser_checkbox, 'Value')
        if (get(handles.blue_shutter, 'Value') == 0)
            set(handles.blue_shutter, 'Value', 1)
        end
    end
    
    if get(handles.red_laser_checkbox, 'Value')
        if (get(handles.red_shutter, 'Value') == 0)
            set(handles.red_shutter, 'Value', 1)
        end
    end
else
    if get(handles.blue_laser_checkbox, 'Value')
        if (get(handles.blue_shutter, 'Value') == 1)
            set(handles.blue_shutter, 'Value', 0)
        end
    end
    
    if get(handles.red_laser_checkbox, 'Value')
        if (get(handles.red_shutter, 'Value') == 1)
            set(handles.red_shutter, 'Value', 0)
        end
    end
end

handles.paramset = 1; %updates parameters
solispreview_Callback(hObject, eventdata, handles)

%turn on/off lasers
digitalLines = [0 get(handles.blue_shutter, 'Value') get(handles.red_shutter, 'Value')]; % Set blue channel to high if button pressed
handles.dioSession.outputSingleScan(digitalLines);

function data_save_path_Callback(hObject, eventdata, handles)
% hObject    handle to data_save_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_save_path as text
%        str2double(get(hObject,'String')) returns contents of data_save_path as a double



function scan_name_suffix_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);


% --- Executes during object creation, after setting all properties.
function scan_name_suffix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scan_name_suffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in galvoright.
function galvoright_Callback(hObject, eventdata, handles)
% hObject    handle to galvoright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
isHR = strcmp('HR', handles.tgroup.SelectedTab.Title);
if(isHR == 1)
    handles.scanFOV = str2double(get(handles.scan_FOV, 'String'));
    offset = str2double(get(handles.galvo_um_offset, 'String'));
    maxPosition = offset+(handles.scanFOV/2);
else
    handles.scanFOV2 = str2double(get(handles.scan_FOV2, 'String'));
    offset = str2double(get(handles.galvo_um_offset, 'String'));
    maxPosition = offset+(handles.scanFOV2/2);
end

set(handles.galvo_position_display, 'String', maxPosition);

galvo_position_display_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function galvo_position_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to galvo_position_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in togglebutton7.
function togglebutton7_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton7


% --- Executes on button press in max_frame_rate.
function checkbox18_Callback(hObject, eventdata, handles)
% hObject    handle to max_frame_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of max_frame_rate



function edit69_Callback(hObject, eventdata, handles)
% hObject    handle to edit69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit69 as text
%        str2double(get(hObject,'String')) returns contents of edit69 as a double


% --- Executes during object creation, after setting all properties.
function edit69_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton6.
function togglebutton6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton6



function edit68_Callback(hObject, eventdata, handles)
% hObject    handle to edit68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit68 as text
%        str2double(get(hObject,'String')) returns contents of edit68 as a double


% --- Executes during object creation, after setting all properties.
function edit68_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit67_Callback(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit67 as text
%        str2double(get(hObject,'String')) returns contents of edit67 as a double


% --- Executes during object creation, after setting all properties.
function edit67_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit66_Callback(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit66 as text
%        str2double(get(hObject,'String')) returns contents of edit66 as a double


% --- Executes during object creation, after setting all properties.
function edit66_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit64_Callback(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit64 as text
%        str2double(get(hObject,'String')) returns contents of edit64 as a double


% --- Executes during object creation, after setting all properties.
function edit64_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit63_Callback(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit63 as text
%        str2double(get(hObject,'String')) returns contents of edit63 as a double


% --- Executes during object creation, after setting all properties.
function edit63_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit62_Callback(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit62 as text
%        str2double(get(hObject,'String')) returns contents of edit62 as a double


% --- Executes during object creation, after setting all properties.
function edit62_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit61_Callback(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit61 as text
%        str2double(get(hObject,'String')) returns contents of edit61 as a double


% --- Executes during object creation, after setting all properties.
function edit61_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit60_Callback(hObject, eventdata, handles)
% hObject    handle to edit60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit60 as text
%        str2double(get(hObject,'String')) returns contents of edit60 as a double


% --- Executes during object creation, after setting all properties.
function edit60_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit59_Callback(hObject, eventdata, handles)
% hObject    handle to edit59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit59 as text
%        str2double(get(hObject,'String')) returns contents of edit59 as a double


% --- Executes during object creation, after setting all properties.
function edit59_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit70_Callback(hObject, eventdata, handles)
% hObject    handle to edit70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit70 as text
%        str2double(get(hObject,'String')) returns contents of edit70 as a double


% --- Executes during object creation, after setting all properties.
function edit70_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in galvoleft.
function galvoleft_Callback(hObject, eventdata, handles)
% hObject    handle to galvoleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
isHR = strcmp('HR', handles.tgroup.SelectedTab.Title);
if(isHR == 1)
    handles.scanFOV = str2double(get(handles.scan_FOV, 'String'));
    offset = str2double(get(handles.galvo_um_offset, 'String'));
    minPosition = offset-(handles.scanFOV/2);
else
    handles.scanFOV2 = str2double(get(handles.scan_FOV2, 'String'));
    offset = str2double(get(handles.galvo_um_offset, 'String'));
    minPosition = offset-(handles.scanFOV2/2);
end

set(handles.galvo_position_display, 'String', minPosition);

galvo_position_display_Callback(hObject, eventdata, handles);



% --- Executes during object creation, after setting all properties.
function camera_framerate2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camera_framerate2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function x_ROI2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_ROI2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function y_ROI2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_ROI2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function camera_left2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camera_left2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pixel_binning2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixel_binning2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function step_size2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_size2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function num_scan_steps2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_scan_steps2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function scan_length2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scan_length2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function scan_FOV2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scan_FOV2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function volumetric_scan_rate2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to volumetric_scan_rate2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function stage_slider_Callback(hObject, eventdata, handles)
% hObject    handle to stage_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

newpos = get(handles.stage_slider, 'Value');
if newpos<0; newpos = 0; end
if newpos>50; newpos = 50; end
set(handles.stage_pos,'String',num2str(newpos))
handles.stage.SetAbsMovePos(0,single(newpos));
handles.stage.MoveAbsolute(0,0);

% verify stage position once stage stops moving
timeout = 10;
t1 = clock; % current time
while(etime(clock,t1)<timeout)
    % Verify stage position
    pos = handles.stage.GetPosition_Position(0);
    if (pos < newpos+0.0005 && pos > newpos-0.0005)
        disp('stage moved')
        break;
    end
end
set(handles.stage_pos, 'String', num2str(pos));
set(handles.stage_slider, 'Value', pos);


% --- Executes during object creation, after setting all properties.
function stage_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stage_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function stage_pos_Callback(hObject, eventdata, handles)
% hObject    handle to stage_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newpos = str2double(get(handles.stage_pos, 'String'));
if newpos<0; newpos = 0; end
if newpos>50; newpos = 50; end
set(handles.stage_slider,'Value',newpos)
handles.stage.SetAbsMovePos(0,single(newpos));
handles.stage.MoveAbsolute(0,0);

% verify stage position once stage stops moving
timeout = 10;
t1 = clock; % current time
while(etime(clock,t1)<timeout)
    % Verify stage position
    pos = handles.stage.GetPosition_Position(0);
    if (pos < newpos+0.0005 && pos > newpos-0.0005)
        disp('stage moved')
        break;
    end
end
set(handles.stage_pos, 'String', num2str(pos));
set(handles.stage_slider, 'Value', pos);





% --- Executes during object creation, after setting all properties.
function stage_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stage_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stage_dist_Callback(hObject, eventdata, handles)
% hObject    handle to stage_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stage_dist as text
%        str2double(get(hObject,'String')) returns contents of stage_dist as a double


% --- Executes during object creation, after setting all properties.
function stage_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stage_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stage_controls.
function stage_controls_Callback(hObject, eventdata, handles)
% hObject    handle to stage_controls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    set(handles.f, 'Visible', 'on')
catch
    waith = waitbar(0, 'Initializing stage. Please wait. And remember to Home before use!');
    set(handles.motorstage,'Value', 1)
    fpos    = get(0,'DefaultFigurePosition'); % figure default position
    fpos(3) = 500; % figure window size;Width
    fpos(4) = 337; % Height
    
    handles.f = figure('Position', fpos,...
        'Menu','None',...
        'Name','APT GUI','Visible','on');
    
    handles.stage = actxcontrol('MGMOTOR.MGMotorCtrl.1',[0 0 500 337 ], handles.f);
    handles.stage.StartCtrl;
    
    % Set the Serial Number
    SN = 28250429; % put in the serial number of the hardware
    set(handles.stage,'HWSerialNum', SN);
    
    % Indentify the device
    handles.stage.Identify;
    for ii = 1:5
        waitbar(ii/5)
        pause(1);
    end
    % Verify stage position
    pos = handles.stage.GetPosition_Position(0);
    set(handles.stage_pos, 'String', num2str(pos));
    if pos<0; pos = 0; elseif pos>50; pos = 50; end % so slider appears at startup when position may not be calibrated
    set(handles.stage_slider, 'Value', pos);
    stepsize = str2num(get(handles.stage_step, 'String'))*0.001; %convert to mm
    set(handles.stage_slider,'SliderStep', [stepsize/50 , 10*stepsize/50]);

    close(waith)
    disp('stage ready to use')
end


function stage_speed_Callback(hObject, eventdata, handles)
% hObject    handle to stage_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stage_speed as text
%        str2double(get(hObject,'String')) returns contents of stage_speed as a double
newvel = str2num(get(handles.stage_speed, 'String'));
if newvel>100
    newval = 100;
    disp('Use a lower speed!')
    set(handles.stage_speed, 'String',newvel)
elseif newvel<0;
    newvel = 0;
    disp('Only positive speeds allowed!')
    set(handles.stage_speed, 'String',newvel)
end

handles.stage.SetVelParams(0,0,1000,newvel);



% --- Executes during object creation, after setting all properties.
function stage_speed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stage_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stage_step_Callback(hObject, eventdata, handles)
% hObject    handle to stage_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stepsize = str2num(get(handles.stage_step, 'String'))*0.001; %convert to mm
set(handles.stage_slider,'SliderStep', [stepsize/50 , 10*stepsize/50]);


% --- Executes during object creation, after setting all properties.
function stage_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stage_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stage_home.
function stage_home_Callback(hObject, eventdata, handles)
% hObject    handle to stage_home (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.stage.MoveHome(0,0);


timeout = 10;
t1 = clock; % current time
while(etime(clock,t1)<timeout)
    % Verify stage position
    pos = handles.stage.GetPosition_Position(0);
    if (pos < 0.0005 && pos > -0.0005)
        disp('stage homed')
        break;
    end
end
set(handles.stage_pos, 'String', num2str(pos));
set(handles.stage_slider, 'Value', pos);



% --- Executes on button press in motorstage.
function motorstage_Callback(hObject, eventdata, handles)
% hObject    handle to motorstage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.motorstage,'Value')
    waith = waitbar(0, 'Initializing stage. Please wait. And remember to Home before use!');
    fpos    = get(0,'DefaultFigurePosition'); % figure default position
    fpos(3) = 500; % figure window size;Width
    fpos(4) = 337; % Height
    
    handles.f = figure('Position', fpos,...
        'Menu','None',...
        'Name','APT GUI','Visible','off');
    
    handles.stage = actxcontrol('MGMOTOR.MGMotorCtrl.1',[0 0 500 337], handles.f);
    handles.stage.StartCtrl;
    
    % Set the Serial Number
    SN = 28250429; % put in the serial number of the hardware
    set(handles.stage,'HWSerialNum', SN);
    
    % Indentify the device
    handles.stage.Identify;
    for ii = 1:5
        waitbar(ii/5)
        pause(1);
    end
    % Verify stage position
    pos = handles.stage.GetPosition_Position(0);
    set(handles.stage_pos, 'String', num2str(pos));
    if pos<0; pos = 0; elseif pos>50; pos = 50; end % so slider appears at startup when position may not be calibrated
    set(handles.stage_slider, 'Value', pos);
    stepsize = str2num(get(handles.stage_step, 'String'))*0.001; %convert to mm
    set(handles.stage_slider,'SliderStep', [stepsize/50 , 10*stepsize/50]);
    close (waith)
    disp('stage ready to use')
    
else
    handles.stage.StopCtrl() 
    pause(1)
    clear handles.f fpos handles.stage
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in solispreview.
function solispreview_Callback(hObject, eventdata, handles)
% hObject    handle to solispreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~get(handles.solispreview,'Value')    
    %shut off red and blue lasers if on
    blueShutterValue = get(handles.blue_shutter, 'Value');
    redShutterValue = get(handles.red_shutter, 'Value');
    if (blueShutterValue == 1)
        set(handles.blue_shutter, 'Value', 0)
    end
    if (redShutterValue == 1)
        set(handles.red_shutter, 'Value', 0)
    end
    digitalLines = [0 0 0]; % Set blue channel to high if button pressed
    handles.dioSession.outputSingleScan(digitalLines);
    system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_stopPreview.ahk');
    system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\returnToGUI.ahk');
    scan_calculations(hObject, handles);
else
    if handles.paramset %only updates SOLIS if parameters have changed
        scan_calculations(hObject, handles);
    else
    system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_startPreview.ahk');
    system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\returnToGUI.ahk');
    end
end


% --- Executes on button press in previewFOV.
function previewFOV_Callback(hObject, eventdata, handles)
disp('not yet implemented')


% --- Executes on button press in lowpower.
function lowpower_Callback(hObject, eventdata, handles)
% hObject    handle to lowpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.lowpower, 'Value')
    set(handles.medpower, 'Value', 0)
    set(handles.highpower, 'Value', 0)
    if get(handles.preview, 'Value')
        set(handles.filter_wheel, 'Value', 4);
        input = 'pos=4';
        fprintf(handles.filterWheelObj, input);
    end
end
% Hint: get(hObject,'Value') returns toggle state of lowpower

% --- Executes on button press in medpower.
function medpower_Callback(hObject, eventdata, handles)
% hObject    handle to medpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.medpower, 'Value')
    set(handles.lowpower, 'Value', 0)
    set(handles.highpower, 'Value', 0)
    if get(handles.preview, 'Value')
        set(handles.filter_wheel, 'Value', 4);
        input = 'pos=4';
        fprintf(handles.filterWheelObj, input);
    end
end
% Hint: get(hObject,'Value') returns toggle state of medpower


% --- Executes on button press in highpower.
function highpower_Callback(hObject, eventdata, handles)
% hObject    handle to highpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.highpower, 'Value')
    set(handles.lowpower, 'Value', 0)
    set(handles.medpower, 'Value', 0)
    if get(handles.preview, 'Value')
        set(handles.filter_wheel, 'Value', 3);
        input = 'pos=3';
        fprintf(handles.filterWheelObj, input);
    end
end
% Hint: get(hObject,'Value') returns toggle state of highpower


% --- Executes on button press in selectROI.
function selectROI_Callback(hObject, eventdata, handles)
% hObject    handle to selectROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ***this is basically the scan_calculations function***
info.dataDirectory = get(handles.data_directory, 'String');
info.scanName = 'SelectROI_temp';

% Step 1: Set ROI values on GUI
info.camera.xROI = 2000;
info.camera.yROI = 400;
info.camera.x_left = 21;

% Step 2: Set calibration factors for nikon 10x)
cFact = handles.savedConversionFactors(1);
info.GUIcalFactors.z_umPerPix = cFact.z_umPerPix;     % Does change as a function of objective magnification
info.GUIcalFactors.y_umPerPix = cFact.y_umPerPix;   % Does change as a function of objective magnification
info.GUIcalFactors.xK_umPerVolt = cFact.xK_umPerVolt; % Doesn't change as a function of objective magnification (Volts to galvo/um)
% Step 3: Set initial scan parameters
% Camera Frame Rate
info.camera.framerate = 100;

% Scan Rate (Volumes/sec)
info.daq.scanRate = 0.9782; %empirically determined for FOV and framerate used

% Scan FOV:
info.daq.scanFOV = 1000;

%binfactor
info.camera.binFactor = 1;

set(handles.galvo_um_offset, 'String', '0')
offset_um  = str2num(get(handles.galvo_um_offset, 'String'));

% Number of frames allocated to each volume (add an extra frame for the
% flyback)
info.daq.pixelsPerLine = 100+1;

% Step 4: Set information regarding scan duration
% Total number of scans to acquire and number of seconds to acquire
info.daq.numberOfScans = 1;
info.daq.scanLength = (info.daq.pixelsPerLine/info.camera.framerate);
info.HR = 1;

% Number of frames acquired by camera
info.camera.kineticSeriesLength =  num2str(ceil(info.daq.scanLength*info.camera.framerate)+100);

% Step 5: Set up scan waveform
info.daq.sawtooth = 1;
info.daq.scanAngle = (info.daq.scanFOV/info.GUIcalFactors.xK_umPerVolt);
info.daq.galvoOffset = offset_um/info.GUIcalFactors.xK_umPerVolt;
globalOffset_volts = handles.globalGalvoOffset/info.GUIcalFactors.xK_umPerVolt;
scanPattern = zeros(info.daq.pixelsPerLine, 1);
maxVal = globalOffset_volts+info.daq.galvoOffset+info.daq.scanAngle/2;
minVal = globalOffset_volts+info.daq.galvoOffset-info.daq.scanAngle/2;
scanPattern(1:end-1) = linspace(minVal, maxVal, info.daq.pixelsPerLine-1);
scanPattern(end) = minVal;
scanPattern = repmat(scanPattern, [info.daq.numberOfScans, 1]);
info.GUIcalFactors.x_umPerPix = info.GUIcalFactors.xK_umPerVolt*info.daq.scanAngle/(info.daq.pixelsPerLine-2);
info.daq.scanWaveform = scanPattern';

% Step 7: Read laser parameters
useBlueLaser = get(handles.blue_laser_checkbox, 'Value');
useRedLaser = get(handles.red_laser_checkbox, 'Value');
info.shutterToggle = [useBlueLaser useRedLaser];
%Modifying laser power info to reflect preview mode changes
info.blue_laser_output_power = 2;
info.red_laser_output_power = 2;

% Step 12: update handles
handles.framerateChanged = addlistener(handles.frameRateFile,'Changed', ...
    @(hObject,eventdata)SCAPE_simpleGalvoGUI_v13('frame_rate_changed',hObject,eventdata,handles));

% Step 13: Write the config file for Zyla
presentDirectory = pwd;
cd(handles.GUI_filepath);
fileName = 'cameraConfig.txt';
% Want to read certain lines and change the value of those lines
fid = fopen(fileName, 'w');
dataDirectory = info.dataDirectory;
indices = find(dataDirectory == '\');

for i = 1:length(indices)-1
    dataDirectory = [dataDirectory(1:indices(i)) '\' dataDirectory(indices(i)+1:end)];
    indices = indices+1;
end
dataDirectory = [dataDirectory '\'];

% % % saves the current scan's time stamp to the info file
info.scanStartTimeApprox = datestr(now);
fprintf(fid, [num2str(info.camera.framerate) '\n']);
fprintf(fid, [num2str(info.camera.xROI) '\n']);
fprintf(fid, [num2str(info.camera.yROI) '\n']);
fprintf(fid, [num2str(info.camera.x_left) '\n']);
fprintf(fid, [num2str(info.camera.binFactor) '\n']);
fprintf(fid, [num2str(info.daq.scanLength) '\n']);
fprintf(fid, [dataDirectory '\n']);
fprintf(fid, [info.scanName '\n']);
fclose(fid);
cd(presentDirectory);

% update AI logging file information
info.daq.logFilePath = [dataDirectory info.scanName '_stim_data.bin'];
info.scanStatus = 'incomplete';

% **this is basically the scan_callback function**
system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_stopPreview.ahk'); %precautionary step
system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_setParameters.ahk');
pause(2);
system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\solis_startAcquisition.ahk');
pause(1.5)
system('C:\Users\Axel-SCAPE\Documents\MATLAB\SCAPE_GUI\NEW_2018\returnToGUI.ahk');
clc

% reset handles.stopButtonFlag
handles.stopButtonFlag = 0;

% move galvo to start pos
handles.aoSession.outputSingleScan(info.daq.scanWaveform(1)')

% set digital line value to LOW keeping shutters closed
handles.dioSession.outputSingleScan([0 0 0]);

% set scan rate
handles.aoSession.Rate = info.camera.framerate;

%put data on AO object for scanning
handles.aoSession.queueOutputData(info.daq.scanWaveform');

info.scanStartTimeApprox = datestr(now);

info.daq.pixelrate = info.camera.framerate;

% open log file
FID = fopen(info.daq.logFilePath, 'w', 'n');

% add listener to AI
lh = handles.aiSession.addlistener('DataAvailable', @(src, event)logData(src, event, FID));

% start AO for scan
handles.aoSession.startBackground;

% start AI
handles.aiSession.startBackground;

% open shutters
handles.dioSession.outputSingleScan([0 info.shutterToggle]);

% send start trigger by setting digital line to HIGH
handles.dioSession.outputSingleScan([1 info.shutterToggle]);

% wait until AO done
handles.aoSession.wait;

% close shutters
handles.dioSession.outputSingleScan([0 0 0]);

% Unclick the buttons
set(handles.blue_shutter, 'Value', 0);
set(handles.red_shutter, 'Value', 0);

% wait until AO done
handles.aoSession.wait;

% stop AI
handles.aiSession.stop;
delete(lh);
fclose(FID);

% update current position display
set(handles.galvo_position_display, 'String', num2str(handles.info.daq.scanWaveform(end) / handles.galvoConversionFactor));
set(handles.galvo_position_slider, 'Value', handles.info.daq.scanWaveform(end) / handles.galvoConversionFactor);

% Center the galvo and turn laser power to zero
% KPedit - added in blue laser
laser_power_r = get(handles.red_laser_output_power, 'String');
laser_power_b = get(handles.blue_laser_output_power, 'String');
set(handles.red_laser_output_power, 'String', '0');
set(handles.blue_laser_output_power, 'String', '0');
center_galvo_Callback(hObject, eventdata, handles)
set(handles.red_laser_output_power, 'String', laser_power_r);
set(handles.blue_laser_output_power, 'String', laser_power_b);


% **Read in scan**
% Creates file path of zyla metadata file (.ini)
zylaInfoFilePath = fullfile(dataDirectory, 'SelectROI_temp', 'acquisitionmetadata.ini');
FID = fopen(zylaInfoFilePath, 'r');
zylaMetaData = fread(FID, '*char')';
fclose(FID);
% Reads the number of frames in a spool file
info.zyla.numDepths = str2double(zylaMetaData(12+(strfind(zylaMetaData, 'AOIHeight = ')):(strfind(zylaMetaData, 'AOIWidth = ')-1)));
info.zyla.numLatPix = str2double(zylaMetaData((11+strfind(zylaMetaData, 'AOIWidth = ')):(strfind(zylaMetaData, 'AOIStride = ')-1)));
info.zyla.imageBytes = str2double(zylaMetaData((17+strfind(zylaMetaData, 'ImageSizeBytes = ')):(strfind(zylaMetaData, '[multiimage]')-1)));
info.zyla.numFramesPerSpool = str2double(zylaMetaData((16+strfind(zylaMetaData, 'ImagesPerFile = ')):end));
startIndex = strfind(zylaMetaData, 'ImageSizeBytes') + length('ImageSizeBytes = ');
info.zyla.ImageSize = str2double(zylaMetaData(startIndex : startIndex + 7));
numFramesPerSpool = info.zyla.numFramesPerSpool;
numDepths = info.zyla.numDepths;
numLatPix = info.zyla.numLatPix;

numSecondsToLoad = info.daq.scanLength;

% Scale data set size by spatial bin factor
scanRate = info.daq.scanRate;               % Volumetric Scan Rate (VPS)
numFrames = info.daq.pixelrate*info.daq.scanLength;% length(info.daq.scanWaveform);  % Total number of frames acquired
numScans = info.daq.numberOfScans;          % Total number of volumes acquired
framesPerScan = info.daq.pixelsPerLine;%floor(info.camera.framerate/info.daq.scanRate);        % Number of frames per volume

% numSpoolfilesPerVolume = ceil(handles.info.info.daq.pixelsPerLine/handles.info.zyla.numFramesPerSpool);
numSpoolfilesPerVolume = ceil(framesPerScan/numFramesPerSpool);
numSpoolfilesPerVolume1 = (framesPerScan/numFramesPerSpool);

numActual_spoolfiles = length(dir([dataDirectory,'/',info.scanName,'/*.dat']));
if numSpoolfilesPerVolume1*scanRate*numSecondsToLoad>numActual_spoolfiles
    numSecondsToLoad = floor(numActual_spoolfiles/(numSpoolfilesPerVolume1*scanRate))-1;
    pf = 'pf';
else
    pf = '';
end


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
filesToLoad = namesOut;


% add extra columns for 2 buffer rows
if (2 == info.camera.binFactor)
    numColumns = numDepths + 1;
    numRows = info.zyla.ImageSize / 2 / numColumns;
    if(mod(numRows, 1) ~= 0)
        numColumns = numDepths+2;
        numRows = info.zyla.ImageSize/2/numColumns;
    end
else
    numColumns = numDepths + 2;
    numRows = numLatPix;
end

% Figures out index for first frame in each volume
frame_starts = 1:framesPerScan:framesPerScan*numScans;

volumesToLoad = 1:round(numSecondsToLoad*scanRate);

try
    startFrameInVolumeToLoad = frame_starts(volumesToLoad);
catch
    try
        startFrameInVolumeToLoad = frame_starts(volumesToLoad(1:end-1));
    catch
        keyboard
    end
end

% Creating the index of every single frames.
numVolumeToLoad = length(startFrameInVolumeToLoad);
endFrameInVolumeToLoad = startFrameInVolumeToLoad+framesPerScan-3;
index = zeros(framesPerScan-2,numVolumeToLoad);

for i = 1:numVolumeToLoad
    index(:,i) = startFrameInVolumeToLoad(i):1:endFrameInVolumeToLoad(i);
end
index = reshape(index,[1,((framesPerScan-2)*numVolumeToLoad)]);

SpoolToLoad = unique(ceil(index/numFramesPerSpool));
numFileToLoad = length(SpoolToLoad);
numFrameToLoad = numVolumeToLoad*framesPerScan;

numPixelsToReshape = numRows * numColumns * numFramesPerSpool;

startSpoolInVolumeToLoad = ceil(startFrameInVolumeToLoad/numFramesPerSpool);
startFrameInSpoolToLoad = mod(startFrameInVolumeToLoad,numFramesPerSpool);

% Load all data into this variable (lateral, depth, scan, time)
clear SCAPE_data
SCAPE_data = zeros(numRows, numColumns, numFramesPerSpool,  numFileToLoad, 'uint16');

frameLoadInSpool = zeros(numFramesPerSpool,numFileToLoad);
for i = 1:numFileToLoad
    frameLoadInSpool(:,i) = ((SpoolToLoad(i)-1)*numFramesPerSpool+1):((SpoolToLoad(i)-1)*numFramesPerSpool+numFramesPerSpool);
end
frameLoadInSpool = reshape(frameLoadInSpool,[1,numFramesPerSpool*numFileToLoad]);
index = find(ismember(frameLoadInSpool,index));


for spoolFileCounter = 1:numFileToLoad
    loadfile = filesToLoad{SpoolToLoad(spoolFileCounter)};
    filePath = fullfile(dataDirectory, info.scanName, loadfile);
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
    
    SCAPE_data(:,:,:,spoolFileCounter) = rawData(1:numRows,1:numColumns,:);
    
end

SCAPE_data = reshape(SCAPE_data,[numRows, numColumns, numFramesPerSpool*numFileToLoad]);
SCAPE_data = reshape(SCAPE_data(:,:,index),[numRows, numColumns, framesPerScan-2, numVolumeToLoad]);

ss = size(SCAPE_data);
%y = ([1:ss(1)]-1)*info.GUIcalFactors.y_umPerPix;
y = ([1:ss(1)]-1);
try
    x = ([1:ss(3)]-1)*info.GUIcalFactors.x_umPerPix;
catch
    info.GUIcalFactors.x_umPerPix = info.GUIcalFactors.xK_umPerVolt*info.daq.scanAngle/(handles.info.daq.pixelsPerLine-2);
    x = ([1:ss(3)]-1)*info.GUIcalFactors.x_umPerPix;
end
%z = ([1:ss(2)]-1)*info.GUIcalFactors.z_umPerPix;
z = ([1:ss(2)]-1);

figure(102)
imagesc(y,x,log(double(squeeze(max(SCAPE_data(:,5:end-5,:),[],2))')));
axis image
xlabel('Y (pixels)')
ylabel('X (um)')
title('X-Y MIP (LOG) - select rectangular ROI and double-click in center when finished')
colormap gray
h = imrect(gca, [50 50 numRows-100 (framesPerScan*info.GUIcalFactors.x_umPerPix)-100]);
xypos = wait(h);
delete(h);
left = floor(xypos(1))+info.camera.x_left;
width = ceil(xypos(3));
scanFOV = round(xypos(4));
galvo_offset = abs((scanFOV-info.daq.scanFOV)/2)-xypos(2);


imagesc(y,z,log(double(squeeze(max(SCAPE_data(:,:,:),[],3))')));
axis image
colormap gray
xlabel('Y (pixels)')
ylabel('Z (pixels)')
title('Y-Z MIP (LOG) - select rectangular ROI and double-click in center when finished')
colormap gray
h = imrect(gca, [left-info.camera.x_left 50 width numColumns-100]);
xypos = wait(h);
delete(h);
left = floor(xypos(1));
width = ceil(xypos(3));

% adjust height to center around center of camera chip with max required vertical offset
height = (info.camera.yROI-xypos(4))/2-xypos(2)+xypos(4);


isHR = strcmp('HR', handles.tgroup.SelectedTab.Title);
if(isHR == 1)
    %HR scan parameters
    set(handles.scan_FOV, 'String', num2str(scanFOV));
    set(handles.x_ROI, 'String', num2str(width));
    set(handles.y_ROI, 'String', num2str(height));
    set(handles.camera_left, 'String', num2str(left));
    scan_FOV_Callback(hObject, eventdata, handles)
        
else
    %fast scan parameters
    set(handles.scan_FOV2, 'String', num2str(scanFOV));
    set(handles.x_ROI2, 'String', num2str(width));
    set(handles.y_ROI2, 'String', num2str(height));
    set(handles.camera_left2, 'String', num2str(left));
    scan_FOV2_Callback(hObject, eventdata, handles)
end

set(handles.galvo_um_offset, 'String', num2str(galvo_offset));

close
%delete temp file
rmdir([dataDirectory '\SelectROI_temp'], 's');

% Update handles structure
guidata(hObject, handles);

scan_calculations(hObject,handles);
center_galvo_Callback(hObject, eventdata, handles)


% --- Executes on selection change in alignpower.
function alignpower_Callback(hObject, eventdata, handles)
% hObject    handle to alignpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.preview, 'Value')
userChoice = num2str(get(handles.alignpower, 'Value'));
input = ['pos=' userChoice];
fprintf(handles.filterWheelObj, input);
end
  

% --- Executes during object creation, after setting all properties.
function alignpower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alignpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
