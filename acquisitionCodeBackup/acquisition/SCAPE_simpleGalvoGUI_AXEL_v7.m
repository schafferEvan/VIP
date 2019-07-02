function varargout = SCAPE_simpleGalvoGUI_AXEL_v7(varargin)
% SCAPE_SIMPLEGALVOGUI_AXEL_V7 MATLAB code for SCAPE_simpleGalvoGUI_AXEL_v7.fig
%      SCAPE_SIMPLEGALVOGUI_AXEL_V7, by itself, creates a new SCAPE_SIMPLEGALVOGUI_AXEL_V7 or raises the existing
%      singleton*.
%
%      H = SCAP194/1.57E2_GALVO_GUI_V15 returns the handle to a new SCAPE_SIMPLEGALVOGUI_AXEL_V7 or the handle to
%      the existing singleton*.
%
%      SCAPE_SIMPLEGALVOGUI_AXEL_V7('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCAPE_SIMPLEGALVOGUI_AXEL_V7.M with the given input arguments.
%
%      SCAPE_SIMPLEGALVOGUI_AXEL_V7('Property','Value',...) creates a new SCAPE_SIMPLEGALVOGUI_AXEL_V7 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SCAPE_simpleGalvoGUI_AXEL_v7_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SCAPE_simpleGalvoGUI_AXEL_v7_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDAT4 GUIHANDLES

% Edit the above text to modify the response to help SCAPE_simpleGalvoGUI_AXEL_v7


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SCAPE_simpleGalvoGUI_AXEL_v7_OpeningFcn, ...
                   'gui_OutputFcn',  @SCAPE_simpleGalvoGUI_AXEL_v7_OutputFcn, ...
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

% --- Executes just before SCAPE_simpleGalvoGUI_AXEL_v7 is made visible.
function SCAPE_simpleGalvoGUI_AXEL_v7_OpeningFcn(hObject, eventdata, handles, varargin)
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

% create data acquisition session for dio
handles.dioSession = daq.createSession('ni');
% add digital lines to session

handles.dioSession.addDigitalChannel('Dev1', 'Port0/Line0', 'OutputOnly');
handles.dioSession.addDigitalChannel('Dev1', 'Port0/Line2:4', 'OutputOnly');


% create data acquisition session for AO for scanning
handles.aoSession = daq.createSession('ni');
% create AO object for scanning adding analog output channel 0
handles.aoSession.addAnalogOutputChannel('Dev1', [0:2], 'Voltage');
% set AO for scanning to external triggering
handles.aoSession.addTriggerConnection('external', 'Dev1/PFI0', 'StartTrigger');
handles.aoSession.Connections(1).TriggerCondition = 'RisingEdge';
% set AO for scanning to external clock source on PFI9
handles.aoSession.addClockConnection('external', 'Dev1/PFI9', 'ScanClock');

% create data acquisition session for external exposure signal
handles.aoExposureSession = daq.createSession('ni');
% Add a counter channel to Channel ID 0 aka ctr0. Forces that channel to
% output an 
handles.aoExposureSession.addCounterOutputChannel('Dev1', 0, 'PulseGeneration');
% Set external exposure signal to external triggering
%handles.aoSession.addTriggerConnection('external', 'Dev1/PFI0', 'StartTrigger');
handles.aoSession.Connections(1).TriggerCondition = 'RisingEdge';

% create data acquisition session for AI
handles.aiSession = daq.createSession('ni');
% create AI object adding analog input channels 0:3
aiChannels = str2num( get(handles.ai_input_channels, 'String'));
handles.aiSession.addAnalogInputChannel('Dev1', [0:aiChannels-1], 'Voltage');
% set AI to external triggering
handles.aiSession.addTriggerConnection('external', 'Dev1/PFI0', 'StartTrigger');
handles.aiSession.Connections(1).TriggerCondition = 'RisingEdge';
% 
% % create data acquisition session for AO for stim
% handles.aoStimSession = daq.createSession('ni');
% % create AO object for stimming adding analog output channel 0
% handles.aoStimSession.addAnalogOutputChannel('Dev1', 'ao0', 'Voltage');
% % set AO for stimming to external triggering
% handles.aoStimSession.addTriggerConnection('external', 'Dev1/PFI0', 'StartTrigger');
% % set AO for outputting to external clock source on PFI9
% handles.aoStimSession.addClockConnection('external', 'Dev1/PFI9', 'ScanClock');

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
catch
    obj = instrfind('Name', 'Serial-COM6');
    delete(obj)
    handles.filterWheelObj = serial('COM6', 'BaudRate', 115200, 'DataBits', 8, ...
        'Terminator', 'CR');
    fopen(handles.filterWheelObj);
end


% Get calibration factors
numObjectives = size(get(handles.objective_list, 'String'),1);
handles.calibration_filepath = [handles.GUI_filepath 'SCAPE_calibrations.txt'];
if exist(handles.calibration_filepath,'file')==2
    fid = fopen(handles.calibration_filepath,'r');
    a = fgetl(fid); % Skip 1st line
    a = fgetl(fid); % Skip 2nd line
    
    handles.savedConversionFactors = struct([]);
    for i = 1:numObjectives
        calFactors = fgetl(fid);
        p = strfind(calFactors,'	');    
        handles.savedConversionFactors(i).y_umPerPix = str2num(calFactors(p(2)+1:p(3)-1));
        handles.savedConversionFactors(i).z_umPerPix = str2num(calFactors(p(3)+1:p(4)-1));
        handles.savedConversionFactors(i).xK_umPerVolt = str2num(calFactors(p(4)+1:end));
    end
else
    errordlg('Cannot find SCAPE Calibration file');
end
fclose(fid);

% I change this parameter when I initially set up the system.
handles.globalGalvoOffset = 0;

% Set up listener for when the frame rate update config file changes
handles.frameRateFile = System.IO.FileSystemWatcher(handles.GUI_filepath);
handles.frameRateFile.Filter = 'cameraConfig_frameRateUpdate.txt';
handles.frameRateFile.EnableRaisingEvents = true;

% Set initial parameter variables
handles.cameraFramerate = str2num(get(handles.camera_framerate, 'String'));
handles.framerate_newInput = 1;

handles.scanRate = str2num(get(handles.volumetric_scan_rate, 'String'));
handles.scanFOV = str2num(get(handles.scan_FOV, 'String'));
handles.numScanSteps = str2num(get(handles.num_scan_steps, 'String'));
handles.stepSize = str2num(get(handles.step_size, 'String'));
handles.scanLength = str2num(get(handles.scan_length, 'String'));
% Choose default command line output for SCAPE_simpleGalvoGUI_AXEL_v7
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% find data directory
browse_for_data_directory_Callback(hObject, eventdata, handles)

% UIWAIT makes SCAPE_simpleGalvoGUI_AXEL_v7 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = SCAPE_simpleGalvoGUI_AXEL_v7_OutputFcn(hObject, eventdata, handles) 
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
info.scanName = get(handles.scan_name, 'String');

% Step 1: Set ROI values on GUI
isFullFrame = get(handles.fullFrame, 'Value');
if (isFullFrame == 0)
    info.camera.xROI = str2num(get(handles.x_ROI, 'String'));
    info.camera.yROI = str2num(get(handles.y_ROI, 'String'));
    info.camera.x_left = str2num(get(handles.camera_left, 'String'));
    temp1 = get(handles.pixel_binning, 'String');
    stimPattern = get(handles.pixel_binning, 'Value');
    temp3 = cell2mat(temp1(stimPattern));
    info.camera.binFactor = str2num(temp3(1));
    clear temp1 temp2 temp3
else
    info.camera.xROI = 2042;
    info.camera.yROI = 2042;
    info.camera.x_left = 1;
    temp1 = get(handles.pixel_binning, 'String');
    stimPattern = get(handles.pixel_binning, 'Value');
    temp3 = cell2mat(temp1(stimPattern));
    info.camera.binFactor = str2num(temp3(1));
    clear temp1 temp2 temp3    
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
    cFact = handles.savedConversionFactors(2);
    info.GUIcalFactors.z_umPerPix = cFact.z_umPerPix;     % Does change as a function of objective magnification
    info.GUIcalFactors.y_umPerPix = cFact.y_umPerPix;   % Does change as a function of objective magnification
    info.GUIcalFactors.xK_umPerVolt = cFact.xK_umPerVolt; % Doesn't change as a function of objective magnification (Volts to galvo/um)
end

% Step 3: Set initial scan parameters
% Camera Frame Rate

if (handles.framerate_newInput == 1)
    disp('New fps input!!!')
    info.loadParameters.desiredCameraFramerate = handles.cameraFramerate;
    handles.framerate_newInput = 0;
end
handles.cameraFramerate = str2num(get(handles.camera_framerate, 'String'));
info.camera.framerate = handles.cameraFramerate;
% info.camera.framerate = str2num(get(handles.camera_framerate, 'String'));
info.loadParameters.maxFrameRate = get(handles.max_frame_rate, 'Value');

% Scan Rate (Volumes/sec)
info.daq.scanRate = handles.scanRate;
% Scan FOV:
info.daq.scanFOV = handles.scanFOV;
offset_um  = str2num(get(handles.galvo_um_offset, 'String'));


% Number of frames allocated to each volume (add an extra frame for the
% flyback)
info.daq.pixelsPerLine = handles.numScanSteps+1;


% Step 4: Set information regarding scan duration
% Total number of scans to acquire and number of seconds to acquire
info.loadParameters.isHR = get(handles.HR_checkbox, 'Value');
if (info.loadParameters.isHR == 1)
    info.daq.numberOfScans = 1;
    info.daq.scanLength = (info.daq.pixelsPerLine/info.camera.framerate);
    
    set(handles.scan_length, 'Enable', 'off');
    set(handles.scan_length, 'BackgroundColor', [0.5 0.5 0.5]);
else
    timePerVolume = (info.daq.pixelsPerLine/info.camera.framerate);
    scanLength = handles.scanLength;
    info.daq.numberOfScans = floor(scanLength/timePerVolume);
    if ((scanLength/timePerVolume)<1)
        disp('Must take at least 1 volume');
        info.daq.numberOfScans = 1;
    end
    handles.scanLength  =  info.daq.numberOfScans*timePerVolume;
    info.daq.scanLength = handles.scanLength;
    info.daq.numberOfScans*timePerVolume;
    
    set(handles.scan_length, 'Enable', 'on');
    set(handles.scan_length, 'BackgroundColor', [1 1 1]);
    set(handles.scan_length, 'String', num2str(info.daq.scanLength));
end

% Number of frames acquired by camera
info.camera.kineticSeriesLength =  num2str(ceil(info.daq.scanLength*info.camera.framerate)+100);
info.camera.kineticSeriesLength
% keyboard


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

% Step 6: Set up stim waveform
% keyboard
info.loadParameters.useStim = get(handles.use_stim_checkbox, 'Value');
t = ([1:length(scanPattern)]-1)/info.camera.framerate;
diffTime = mean(diff(t));

if ((info.loadParameters.useStim==1)&&(info.loadParameters.isHR==0))
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
    stimPattern = stimPattern(1:length(scanPattern));
    plot(handles.stim_axes, t, stimPattern);
    xlim([t(1) t(end)]);
else
    stimPattern = scanPattern*0;
    info.daq.stim.Vout = 0;
    cla(handles.stim_axes);
end


% Step 7: Read laser parameters
info.shutterToggle = [1 0 0];
info.laser_output_power = str2num(get(handles.laser_output_power, 'String'));
voltPower = info.laser_output_power*0.0334-0.004;
laserPattern = stimPattern;
laserPattern(2:end) = laserPattern(1:end-1)*0+voltPower;
laserPattern(1) = 0;
info.daq.scanWaveform = [stimPattern'; scanPattern'; laserPattern'];


% Step 8: Update AI channels, sample rate and scan duration
info.daq.aiChannels = [0:str2num( get(handles.ai_input_channels, 'String'))-1];
while 1 <= length(handles.aiSession.Channels)
    handles.aiSession.removeChannel(1);
end
handles.aiSession.addAnalogInputChannel('Dev1', info.daq.aiChannels, 'Voltage');
info.daq.aiSampleRate = str2double( get(handles.ai_sample_rate, 'String'));
handles.aiSession.Rate = info.daq.aiSampleRate;
handles.aiSession.DurationInSeconds = 1 + info.daq.scanLength;

% Step 9:
set(handles.data_save_path, 'String', [info.dataDirectory info.scanName]);

% Step 10: Update the laser power
temp1 = get(handles.filter_wheel, 'String');
stimPattern = get(handles.filter_wheel, 'Value');
info.laser_power = cell2mat(temp1(stimPattern));

% Step 11: Update experiment notes
info.experiment_notes = get(handles.experiment_notes, 'String');

% Step 12: update handles
handles.info = info;
handles.framerateChanged = addlistener(handles.frameRateFile,'Changed', ...
    @(hObject,eventdata)SCAPE_simpleGalvoGUI_AXEL_v7('frame_rate_changed',hObject,eventdata,handles));

% Step 13: Write the config file for Zyla
writeConfigFile(hObject, handles);

guidata(hObject, handles);

disp('Scan Calculations Complete');

function scan_Callback(hObject, eventdata, handles)

cla(handles.stim_axes);
scan_calculations(hObject, handles);

% reset handles.stopButtonFlag
handles.stopButtonFlag = 0;

% move galvo to start pos

handles.aoSession.outputSingleScan(handles.info.daq.scanWaveform(:, 1)');

% 
% % set stim line to 0
% if (1 == get(handles.use_stim_checkbox, 'Value'))
%     handles.aoStimSession.outputSingleScan(handles.info.daq.stim.Vout(1));
% end

% set digital line value to LOW keeping shutters closed
handles.dioSession.outputSingleScan([0 0 0 0]);

% set scan rate
handles.aoSession.Rate = handles.info.camera.framerate;
% 
% set stim output rate
% handles.aoStimSession.Rate = handles.info.daq.stim.fout;

%put data on AO object for scanning
handles.aoSession.queueOutputData(handles.info.daq.scanWaveform');
% 
% % put data on AO object for stimming
% if (1 == get(handles.use_stim_checkbox, 'Value'))
%     handles.aoStimSession.queueOutputData(handles.info.daq.stim.Vout');
% end

% save info file
info = handles.info;
info.camera.kineticSeriesLength = num2str(info.camera.framerate*info.daq.scanLength+100);
info.camera.framerate = str2num((get(handles.camera_framerate, 'String')));
info.daq.pixelrate = info.camera.framerate;
infoFilePath = [info.dataDirectory info.scanName '_info.mat'];
% checks for repeat file names to prevent overwrite of previous info files
fileCounter = 1;
while (2 == exist(infoFilePath, 'file'))
    infoFilePath = [info.dataDirectory info.scanName '_' num2str(fileCounter) '_info.mat'];
    fileCounter = fileCounter + 1;
end

% saves the current scan's time stamp to the info file
info.scanStartTimeApprox = datestr(now);
if (1 < fileCounter)
    info.scanName = [info.scanName '_' num2str(fileCounter - 1)];
end

% update AI logging file information
info.daq.logFilePath = [info.dataDirectory info.scanName '_stim_data.bin'];
info.scanStatus = 'incomplete';

% save preliminary info file
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
% 
% % start AO for stimming
% if (1 == get(handles.use_stim_checkbox, 'Value'))
%     handles.aoStimSession.startBackground;
% end

% start AI
handles.aiSession.startBackground;

% open shutters
handles.dioSession.outputSingleScan([0 1 1 1]);
handles.info.shutterToggle;

% send start trigger by setting digital line to HIGH
tic;
handles.dioSession.outputSingleScan([1 1 1 1]);

% wait until AO done
handles.aoSession.wait;

scanTime = toc;

% close shutters
handles.dioSession.outputSingleScan([0 0 0 0]);
% Unclick the buttons
set(handles.blue_shutter, 'Value', 0)
% set(handles.green_shutter, 'Value', 0)
% set(handles.red_shutter, 'Value', 0)

% wait until AO stim session done
if (1 == get(handles.use_stim_checkbox, 'Value'))
%     handles.aoStimSession.wait;
end

disp(['Done scan! Took: ' num2str(scanTime) ' seconds']);

% wait until AO done
handles.aoSession.wait;

% stop AI
handles.aiSession.stop;
delete(lh);
fclose(FID);

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
laser_power = get(handles.laser_output_power, 'String');
set(handles.laser_output_power, 'String', '0');
center_galvo_Callback(hObject, eventdata, handles)
set(handles.laser_output_power, 'String', laser_power);


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

%
handles = rmfield(handles,'lastStimData'); % Beth added 12/18/14

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code + scan calculations callback %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fullFrame_Callback(hObject, eventdata, handles)
fullFrame = get(handles.fullFrame, 'Value');
if (fullFrame == 1)
   set(handles.x_ROI, 'Enable', 'off');
   set(handles.x_ROI, 'BackgroundColor', [0.5 0.5 0.5]);
   set(handles.y_ROI, 'Enable', 'off');
   set(handles.y_ROI, 'BackgroundColor', [0.5 0.5 0.5]);
   set(handles.camera_left, 'Enable', 'off');
   set(handles.camera_left, 'BackgroundColor', [0.5 0.5 0.5]);
else
    set(handles.x_ROI, 'Enable', 'on');
    set(handles.x_ROI, 'BackgroundColor', [1 1 1]);
    set(handles.y_ROI, 'Enable', 'on');
    set(handles.y_ROI, 'BackgroundColor', [1 1 1]);
    set(handles.camera_left, 'Enable', 'on');
    set(handles.camera_left, 'BackgroundColor', [1 1 1]);
end
scan_calculations(hObject, handles);

function load_info_file_Callback(hObject, eventdata, handles)
presentDirectory = pwd;
cd('D:\');
[fileName, pathName, filterIndex] = uigetfile('*.mat');
cd (pathName)
load(fileName)
cd(presentDirectory)


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
set(handles.max_frame_rate, 'Value', info.loadParameters.maxFrameRate);

% Set volume rate
set(handles.volumetric_scan_rate, 'String', info.daq.scanRate);
handles.scanRate = info.daq.scanRate;

% Set Run Length:
numVolumes = info.daq.numberOfScans;
timePerVolume = 1/handles.scanRate;
handles.scanLength = numVolumes/handles.scanRate +(timePerVolume*0.1);
set(handles.scan_length, 'String', num2str(handles.scanLength));
set(handles.HR_checkbox, 'Value', info.loadParameters.isHR);

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

% Laser power
set(handles.laser_output_power, 'String', num2str(info.laser_output_power));

% Filter wheel
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
fprintf(handles.filterWheelObj, input);
scan_calculations(hObject, handles);

function browse_for_data_directory_Callback(hObject, eventdata, handles)
handles.info.directory = [uigetdir('D:\','Select Directory to save data to') '\'];
cd(handles.info.directory);
% For some reason, any directory chosen that is not a drive (C:, D:, etc)
% is not returned with a backslash, but if a drive is chosen, a backslash
% is appended.
if (handles.info.directory(end-1:end) == '\\')
    handles.info.directory = handles.info.directory(1:end-1)
end

if 0 == handles.info.directory
    errordlg('You did not select a directory. Please choose again.','yeah')
    return;
end

set(handles.data_directory, 'String', handles.info.directory);
guidata(hObject,handles);

scan_calculations(hObject, handles);

function laser_output_power_Callback(hObject, eventdata, handles)
laserPower = str2num(get(handles.laser_output_power, 'String'))/1000;


laserPower = str2num(get(handles.laser_output_power, 'String'));
voltPower = laserPower*0.0334-0.004;

shutterValue = get(handles.blue_shutter, 'Value');
 % get new position
conversionFactor = handles.info.GUIcalFactors.xK_umPerVolt;
newPosition = get(handles.galvo_position_slider, 'Value');
globalOffset = handles.globalGalvoOffset; % Value in microns (set by me when I install the system)

if (1 == shutterValue)
    outputVoltage = [0 (newPosition+globalOffset)/conversionFactor voltPower];
    disp('Blue Shutter Open');
elseif (0 == shutterValue)
    outputVoltage = [0 (newPosition+globalOffset)/conversionFactor 0];
    
    disp('Blue Shutter Closed');
end
handles.aoSession.outputSingleScan(outputVoltage);

scan_calculations(hObject, handles);

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
handles.scanRate = framerate/numFramesPerVolume
set(handles.volumetric_scan_rate, 'String', num2str(handles.scanRate));

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

scan_calculations(hObject, handles);

function scan_length_Callback(hObject, eventdata, handles)
handles.scanLength = str2num(get(handles.scan_length, 'String'));
scan_calculations(hObject, handles);

function frame_rate_changed(hObject, eventdata, handles)
% Update camera frame rate
fileName = [handles.GUI_filepath, 'cameraConfig_frameRateUpdate.txt'];
fid = fopen(fileName);
handles.cameraFramerate = str2num(strtrim(fgets(fid)));
set(handles.camera_framerate, 'String', num2str(handles.cameraFramerate));
fclose(fid);

numScanSteps            = str2num(get(handles.num_scan_steps, 'String'));
% Add frame for flyback
numFramesPerVolume = numScanSteps+1;
% Update parameters
handles.scanRate = handles.cameraFramerate/numFramesPerVolume;
set(handles.volumetric_scan_rate, 'String', num2str(handles.scanRate));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Galvo Related Functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scan_FOV_Callback(hObject, eventdata, handles)
handles.scanFOV = str2double(get(handles.scan_FOV, 'String'));
offset = str2double(get(handles.galvo_um_offset, 'String'));
minPosition = offset-(handles.scanFOV/2);
maxPosition = offset+(handles.scanFOV/2);

% Update step size
info = handles.info;
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
    laserPower = (get(handles.laser_output_power, 'String'));
    set(handles.laser_output_power, 'String', 0);
    galvo_position_display_Callback(hObject, eventdata, handles);
    set(handles.laser_output_power, 'String', laserPower);
else
    galvo_position_display_Callback(hObject, eventdata, handles);
end
    
scan_calculations(hObject, handles);

function galvo_position_display_Callback(hObject, eventdata, handles)

% Get new position
newPosition = str2double( get(handles.galvo_position_display, 'String'));
maxPosition = get(handles.galvo_position_slider, 'Max');
minPosition = get(handles.galvo_position_slider, 'Min');
offset = str2double(get(handles.galvo_um_offset, 'String'));
if ((newPosition > maxPosition) || (newPosition < minPosition))
    % Set new position to max or min value
    if (newPosition > offset)
        newPosition = maxPosition;
    else
        newPosition = -minPosition;
    end
    
    % Update position display
    set(handles.galvo_position_display, 'String', num2str(newPosition, 3));
end

% Update slider with new position value
set(handles.galvo_position_slider, 'Value', newPosition);

% Run slider callback to reposition galvo
galvo_position_slider_Callback(hObject, eventdata, handles);

function center_galvo_Callback(hObject, eventdata, handles)

offset = get(handles.galvo_um_offset, 'String');
set(handles.galvo_position_display, 'String', offset);

galvo_position_display_Callback(hObject, eventdata, handles);

function galvo_um_offset_Callback(hObject, eventdata, handles)
handles.scanFOV = str2double(get(handles.scan_FOV, 'String'));

offset = str2double(get(handles.galvo_um_offset, 'String'));
minPosition = offset-(handles.scanFOV/2);
maxPosition = offset+(handles.scanFOV/2);

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
h = waitbar(0, 'Moving to new position please wait. This will take 1 second');

laserPower = str2num(get(handles.laser_output_power, 'String'));
voltPower = laserPower*0.0334-0.004;

% get new position
conversionFactor = handles.info.GUIcalFactors.xK_umPerVolt;
newPosition = get(handles.galvo_position_slider, 'Value');
globalOffset = handles.globalGalvoOffset; % Value in microns (set by me when I install the system)

% move to position
outputVoltage = [0 (newPosition+globalOffset)/conversionFactor voltPower];
handles.aoSession.outputSingleScan(outputVoltage);

% update position display
handles.currentGalvoPosition = newPosition;
guidata(hObject, handles);
set(handles.galvo_position_display, 'String', num2str(handles.currentGalvoPosition, 3));

% update waitbar
waitbar(1, h);

% close waitbar
close(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Pure Code Functions %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function stop_scan_Callback(hObject, eventdata, handles)
% close shutters
handles.dioSession.outputSingleScan([0 0 0 0]);

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

function ai_channels_to_plot_Callback(hObject, eventdata, handles)
if (1 == isfield(handles, 'lastStimData'))
    channelsToPlot = 1 + str2num( get(handles.ai_channels_to_plot, 'String'));
    plot(handles.stim_axes, handles.lastStimData.time, handles.lastStimData.channels(channelsToPlot, :));
    
    legendText = [];
    if (find(channelsToPlot == 1))
        legendText = {'Galvo Drive'};
    end
    if (find(channelsToPlot == 2))
        legendText = [legendText {'Galvo Feedback'}];
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
for i = 1:length(indices)-1
   dataDirectory = [dataDirectory(1:indices(i)) '\' dataDirectory(indices(i)+1:end)];
   indices = indices+1;
end
dataDirectory = [dataDirectory '\'];


if (handles.info.loadParameters.maxFrameRate == 1)
    framerate = 3000;
else
    framerate = handles.info.loadParameters.desiredCameraFramerate;
end 
    
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
shutterValue = get(handles.blue_shutter, 'Value');
 % get new position
conversionFactor = handles.info.GUIcalFactors.xK_umPerVolt;
newPosition = get(handles.galvo_position_slider, 'Value');
globalOffset = handles.globalGalvoOffset; % Value in microns (set by me when I install the system)
laserPower = str2num(get(handles.laser_output_power, 'String'));
voltPower = laserPower*0.0334-0.004;

if (1 == shutterValue)
    outputVoltage = [0 (newPosition+globalOffset)/conversionFactor voltPower];
    disp('Blue Shutter Open');
elseif (0 == shutterValue)
    outputVoltage = [0 (newPosition+globalOffset)/conversionFactor 0];
    
    disp('Blue Shutter Closed');
end
handles.aoSession.outputSingleScan(outputVoltage);
guidata(hObject, handles);
handles.dioSession.outputSingleScan([0 handles.info.shutterToggle]);

function update_extra_info_parameters(hObject, handles)
% This function saves parameters to an info file that 1) aren't used for
% data loading and 2) are required to load the info file for future
% acquisitions
handles.cameraFramerate = str2num(get(handles.camera_framerate, 'String'));
handles.stepSize = str2num(get(handles.step_size, 'String'));
handles.scanLength = str2num(get(handles.scan_length, 'String'));




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
scan_calculations(hObject, handles);

function scan_name_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function data_directory_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function x_ROI_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function y_ROI_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function camera_binning_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);


function ai_input_channels_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function ai_sample_rate_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function HR_checkbox_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

% --- Executes on button press in max_frame_rate.
function max_frame_rate_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles)

function stimTable_CellEditCallback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function camera_left_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function optical_scan_angle_Callback(hObject, eventdata, handles)
scan_calculations(hObject, handles);

function blue_laser_checkbox_Callback(hObject, eventdata, handles)
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

function red_laser_checkbox_Callback(hObject, eventdata, handles)
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
volumetricScanRate = str2double(get(handles.volumetric_scan_rate, 'String'));
numVolumesAcquired = str2double(get(handles.volumes_to_collect, 'String'));
scanLength = numVolumesAcquired/volumetricScanRate;
set(handles.scan_length, 'String', num2str(scanLength));
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

function red_shutter_Callback(hObject, eventdata, handles)
shutterValue = get(handles.red_shutter, 'Value');
if (1 == shutterValue)
    handles.info.shutterToggle(3) = 1;
    disp('Red Shutter Open');
elseif (0 == shutterValue)
    handles.info.shutterToggle(3) = 0;
    disp('Red Shutter Open');
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
        conversionFactors(2) = cFact.lateral_umPerPix;   % Does change as a function of objective magnification
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
% hObject    handle to laser_output_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
