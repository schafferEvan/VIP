function varargout = SCAPE_quickMovie_ES(varargin)
% SCAPE_QUICKMOVIE_ES MATLAB code for SCAPE_quickMovie_ES.fig
%      SCAPE_QUICKMOVIE_ES, by itself, creates a new SCAPE_QUICKMOVIE_ES or raises the existing
%      singleton*.
%
%      H = SCAPE_QUICKMOVIE_ES returns the handle to a new SCAPE_QUICKMOVIE_ES or the handle to
%      the existing singleton*.
%
%      SCAPE_QUICKMOVIE_ES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCAPE_QUICKMOVIE_ES.M with the given input arguments.
%
%      SCAPE_QUICKMOVIE_ES('Property','Value',...) creates a new SCAPE_QUICKMOVIE_ES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SCAPE_quickMovie_ES_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SCAPE_quickMovie_ES_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SCAPE_quickMovie_ES

% Last Modified by GUIDE v2.5 20-Apr-2018 14:32:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SCAPE_quickMovie_ES_OpeningFcn, ...
                   'gui_OutputFcn',  @SCAPE_quickMovie_ES_OutputFcn, ...
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

% --- Executes just before SCAPE_quickMovie_ES is made visible.
function SCAPE_quickMovie_ES_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SCAPE_quickMovie_ES (see VARARGIN)

% Choose default command line output for SCAPE_quickMovie_ES
handles.output = hObject;
handles.data = varargin{1};
handles.ss = size(handles.data);
handles.info = varargin{2};

% Set default slice numbers
handles.t = 1; 
handles.xslice = 1;
handles.yslice = 1;
handles.zslice = 1;

% Set axes
handles.y = [1:handles.ss(1)]-1;
handles.x = [1:handles.ss(3)]-1;
handles.z = [1:handles.ss(2)]-1;


% Set default colormap
if (length(handles.ss) == 3)
    handles.maxDataColor = double(squeeze(max(max(max(handles.data)))));
elseif (length(handles.ss) == 4)
    handles.maxDataColor = double(squeeze(max(max(max(max(handles.data))))));
end
set(handles.colorMaxTextbox, 'String', num2str(handles.maxDataColor));
% Set default crop parameters
handles.cropParams.xmin = handles.x(1);
handles.cropParams.ymin = handles.y(1);
handles.cropParams.zmin = handles.z(1);
handles.cropParams.xmax = 10000;%handles.x(end);
handles.cropParams.ymax = 10000;%handles.y(end);
handles.cropParams.zmax = 10000; %handles.z(end);
set(handles.cropXmaxTextbox, 'String', num2str(handles.cropParams.xmax))
set(handles.cropYmaxTextbox, 'String', num2str(handles.cropParams.ymax))
set(handles.cropZmaxTextbox, 'String', num2str(handles.cropParams.zmax))

% Update handles structure
axis(handles.axes1);

imagesc(handles.y, handles.x, squeeze(handles.data(:,1,:,1))');
axis image
colormap gray;
xlabel('Y', 'Color', 'w')
ylabel('X', 'Color', 'w')
set(handles.axes1, 'YColor', 'w');
set(handles.axes1, 'XColor', 'w');

if(length(handles.ss) == 3)
    sliderStep = 1;
    set(handles.timeStamp, 'Visible', 'Off');
    set(handles.timeslider, 'Enable', 'Off');
    set(handles.timenum, 'Enable', 'Off');
    set(handles.playPushbutton, 'Enable', 'Off');
    set(handles.slowerPushbutton, 'Enable', 'Off');
    set(handles.fasterPushbutton, 'Enable', 'Off');
    set(handles.timecoursePushbutton, 'Enable', 'Off');
    set(handles.clearTimecoursesPushbutton, 'Enable', 'Off');
    set(handles.exportTimecoursesPushbutton, 'Enable', 'Off');
else
    sliderStep = 1/(handles.ss(4)-1);
end
set(handles.timeslider, 'SliderStep', [sliderStep, sliderStep*5])
sliceSliderStep = 1/(handles.ss(2)-1);
set(handles.sliceslider, 'SliderStep', [sliceSliderStep, sliceSliderStep*5])
handles.playbackSpeed = 100;
handles.numTC = 1;

if isfield(handles.info,'listbox'); set(handles.infobox,'String',handles.info.listbox); end
guidata(hObject, handles);

% UIWAIT makes SCAPE_quickMovie_ES wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = SCAPE_quickMovie_ES_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on slider movement.
function sliceslider_Callback(hObject, eventdata, handles)
% hObject    handle to sliceslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sl = get(handles.sliceslider,'Value');

if (sl==0)
    sl = 0.00000001; % Prevents the slice number from being set to zero.
end

colorMaxVal = get(handles.colorMaxSlider, 'Value')*handles.maxDataColor;
colorMinVal = get(handles.colorMinSlider, 'Value')*handles.maxDataColor;
tmp = uint16([colorMinVal colorMaxVal]);

logscale = get(handles.logScaleCheckbox, 'Value');
if(logscale)
    clims = log10(double([min(tmp) max(tmp)]));
else
    clims = [min(tmp) max(tmp)];
end

mipview = get(handles.mipViewCheckbox, 'Value');
xwin = round(str2num(get(handles.xSliceThicknessTextbox, 'String'))/2);
ywin = round(str2num(get(handles.ySliceThicknessTextbox, 'String'))/2);
zwin = round(str2num(get(handles.zSliceThicknessTextbox, 'String'))/2);
    
if get(handles.xyzview,'Value')==1;  
    if(mipview)
        win = [-zwin zwin];
    else
        win = [0 0];
    end
    sll = ceil(handles.ss(2)*sl);
    handles.zslice = sll+win;
    handles.zslice(1) = max(handles.zslice(1), 1);
    handles.zslice(2) = min(handles.zslice(2), handles.ss(2));
    handles.zslice = [handles.zslice(1):handles.zslice(2)];
    
    axis(handles.axes1);
    if(~logscale)
        imagesc(handles.y, handles.x, squeeze(max(handles.data(:,handles.zslice,:,handles.t), [], 2))');
    else
        imagesc(handles.y, handles.x, log10(double(squeeze(max(handles.data(:,handles.zslice,:,handles.t), [], 2))')));
    end
    axis image
    xlim([max(handles.cropParams.ymin, 0) min(handles.y(end), handles.cropParams.ymax)])
    ylim([max(handles.cropParams.xmin, 0) min(handles.x(end), handles.cropParams.xmax)])
    xlabel('Y (um)', 'Color', 'w')
    ylabel('X (um)', 'Color', 'w')
    
elseif get(handles.xyzview,'Value')==2;
    if(mipview)
        win = [-xwin xwin];
    else
        win = [0 0];
    end
    sll = ceil(handles.ss(3)*sl);
    handles.xslice = sll+win;
    handles.xslice(1) = max(handles.xslice(1), 1);
    handles.xslice(2) = min(handles.xslice(2), handles.ss(3));
    handles.xslice = [handles.xslice(1):handles.xslice(2)];
    
    axis(handles.axes1);
    if(~logscale)
        imagesc(handles.y, handles.z, squeeze(max(handles.data(:,:,handles.xslice,handles.t), [], 3))');
    else
        imagesc(handles.y, handles.z, log10(double(squeeze(max(handles.data(:,:,handles.xslice,handles.t), [], 3))')));
    end
    ylabel('Z (um)', 'Color', 'w')
    xlabel('Y (um)', 'Color', 'w')
    axis image
    xlim([max(handles.cropParams.ymin, 0) min(handles.y(end), handles.cropParams.ymax)])
    ylim([max(handles.cropParams.zmin, 0) min(handles.z(end), handles.cropParams.zmax)])
    
elseif get(handles.xyzview,'Value')==3;    
    if(mipview)
        win = [-ywin ywin];
    else
        win = [0 0];
    end
    sll = ceil(handles.ss(1)*sl);
    handles.yslice = sll+win;
    handles.yslice(1) = max(handles.yslice(1), 1);
    handles.yslice(2) = min(handles.yslice(2), handles.ss(1));
    handles.yslice = [handles.yslice(1):handles.yslice(2)];
    
    axis(handles.axes1);
    if(~logscale)
        imagesc(handles.x, handles.z, squeeze(max(handles.data(handles.yslice,:,:,handles.t), [], 1)));
    else
        imagesc(handles.x, handles.z, log10(double(squeeze(max(handles.data(handles.yslice,:,:,handles.t), [], 1)))));
    end
    ylabel('Z (um)', 'Color', 'w')
    xlabel('X (um)', 'Color', 'w')
    axis image
    xlim([max(handles.cropParams.xmin, 0) min(handles.cropParams.xmax, handles.x(end))])
    ylim([max(handles.cropParams.zmin, 0) min(handles.z(end), handles.cropParams.zmax)])
end
set(handles.timeStamp, 'String', ['Time: ' num2str(round(handles.t/handles.info.daq.scanRate*100)/100) ' sec']);
caxis(clims);
colormap gray;
set(handles.axes1, 'YColor', 'w');
set(handles.axes1, 'XColor', 'w');
set(handles.slicenum,'String',sll);
guidata(hObject, handles);


function slicenum_Callback(hObject, eventdata, handles)
% hObject    handle to slicenum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mipval = get(handles.xyzview,'Value');
sliceval = str2num(get(handles.slicenum, 'String'));

switch mipval
    case 1
        maxsliceval = handles.ss(2);
    case 2
        maxsliceval = handles.ss(3);
    case 3
        maxsliceval = handles.ss(1);
end
if(sliceval> maxsliceval)
    sliceval = maxsliceval;
    set(handles.slicenum, 'String', num2str(sliceval));
end
if (sliceval < 1)
    sliceval = 1;
    set(handles.slicenum, 'String', num2str(sliceval));
end
set(handles.sliceslider,'Value', sliceval/maxsliceval);

sliceslider_Callback(hObject, eventdata, handles)

% --- Executes on selection change in xyzview.
function xyzview_Callback(hObject, eventdata, handles)
% hObject    handle to xyzview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xyzview contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xyzview
val = get(handles.xyzview, 'Value');

switch val
    case 1 
        sliceSliderStep = 1/(handles.ss(2)-1);
        set(handles.sliceslider, 'SliderStep', [sliceSliderStep, sliceSliderStep*5])
        set(handles.xSliceThicknessTextbox, 'Enable','Off');
        set(handles.ySliceThicknessTextbox, 'Enable','Off');
        set(handles.zSliceThicknessTextbox, 'Enable','On');
    case 2
        sliceSliderStep = 1/(handles.ss(3)-1);
        set(handles.sliceslider, 'SliderStep', [sliceSliderStep, sliceSliderStep*5])
        set(handles.xSliceThicknessTextbox, 'Enable','On');
        set(handles.ySliceThicknessTextbox, 'Enable','Off');
        set(handles.zSliceThicknessTextbox, 'Enable','Off');
    case 3
        sliceSliderStep = 1/(handles.ss(1)-1);
        set(handles.sliceslider, 'SliderStep', [sliceSliderStep, sliceSliderStep*5])
        set(handles.xSliceThicknessTextbox, 'Enable','Off');
        set(handles.ySliceThicknessTextbox, 'Enable','On');
        set(handles.zSliceThicknessTextbox, 'Enable','Off');
end


sliceslider_Callback(hObject, eventdata, handles)

%%%
%%% TIME PLAYBACK CALLBACKS
%%%

function timeslider_Callback(hObject, eventdata, handles)
% hObject    handle to timeslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if (length(handles.ss) == 3)
    handles.t = 1;
else
    handles.t = ceil(handles.ss(4)*get(handles.timeslider,'Value'));
    if(handles.t == 0)
        handles.t = 1;
    end
end
set(handles.timenum,'String',num2str(handles.t));
sliceslider_Callback(hObject, eventdata, handles)

function timenum_Callback(hObject, eventdata, handles)
% hObject    handle to timenum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = str2num(get(handles.timenum, 'String'));
if (val <=0)
    val = 1
elseif (val >= handles.ss(4))
    val = handles.ss(4);
end
handles.t = val;
set(handles.timeslider, 'Value', val/handles.ss(4));
timeslider_Callback(hObject, eventdata, handles)


function movHandle = setupMovie(info)
slashes = strfind(info.dataDirectory,'\');
movName = [info.dataDirectory(1:slashes(end-1)),'_movies',info.dataDirectory(slashes(end-1):end-1),'_',info.scanName];
movHandle = VideoWriter(movName);
movHandle.FrameRate = 12;
open(movHandle);


function playPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to playPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(handles.playPushbutton, 'Value');
set(handles.timeslider, 'Enable', 'Off');
set(handles.slowerPushbutton, 'Enable', 'Off');
set(handles.fasterPushbutton, 'Enable', 'Off');
vid = setupMovie(handles.info);

while(val == 1)
    val = get(handles.playPushbutton, 'Value');
    handles.t = handles.t+1;
    if(handles.t>handles.ss(4))
        handles.t = 1;
    end
    set(handles.timenum,'String',num2str(handles.t));
    sliceslider_Callback(hObject, eventdata, handles);
    writeVideo(vid, getframe(handles.axes1));
    drawnow
    pause(max(1/handles.playbackSpeed,.1))
end
close(vid);
set(handles.timeslider, 'Enable', 'On');
set(handles.slowerPushbutton, 'Enable', 'On');
set(handles.fasterPushbutton, 'Enable', 'On');
timenum_Callback(hObject, eventdata, handles)

function slowerPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to slowerPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.playbackSpeed = handles.playbackSpeed*0.8;
guidata(hObject, handles);

function fasterPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to fasterPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.playbackSpeed = handles.playbackSpeed*1.2;
guidata(hObject, handles);


% --- Executes on button press in timecoursePushbutton.
function timecoursePushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to timecoursePushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mipval = get(handles.xyzview, 'Value');
calval = (get(handles.uniformAspectRadiobutton, 'Value')*1)+...
    (get(handles.calibratedAspectRadiobutton, 'Value')*2)+...
    (get(handles.customAspectRadiobutton, 'Value')*3);
slicenum = str2num(get(handles.slicenum, 'String'));
cmap = lines(50);
hold on;
switch mipval
    case 1
        [y,x] = ginput_col(2);
        x = round(x);
        y = round(y);
        z = slicenum;
        plot(handles.axes1, mean(y), mean(x), 'o', 'color', cmap(handles.numTC, :));
    case 2
        [y,z] = ginput_col(2);
        y = round(y);
        z = round(z);
        x = slicenum;
        plot(handles.axes1, mean(y), mean(z), 'o', 'color', cmap(handles.numTC, :));
    case 3
        [x,z] = ginput_col(2);
        x = round(x);
        z = round(z);
        y = slicenum;
        plot(handles.axes1, mean(x), mean(z), 'o', 'color', cmap(handles.numTC, :));
end
hold off

time = ([1:handles.ss(4)]-1)/handles.info.daq.scanRate;
switch calval
    case 1
        % Do nothing. The axes coordinates and the indices of the data
        % matrix are already the same.
    case 2
        x = (x/handles.info.GUIcalFactors.x_umPerPix)+1;
        y = (y/handles.info.GUIcalFactors.y_umPerPix)+1;
        z = (z/handles.info.GUIcalFactors.z_umPerPix)+1;
    case 3
        xCal = str2num(get(handles.xCalTextbox, 'String'));
        yCal = str2num(get(handles.yCalTextbox, 'String'));
        zCal = str2num(get(handles.zCalTextbox, 'String'));
        x = (x/xCal)+1;
        y = (y/yCal)+1;
        z = (z/zCal)+1;
end

switch mipval
    case 1
        z = slicenum;
    case 2
        x = slicenum;
    case 3
        y = slicenum;
end


x =round(x); y = round(y); z = round(z);
handles.ROI(handles.numTC).TC = squeeze(mean(mean(squeeze(handles.data(min(y):max(y), ...
    min(z):max(z), min(x):max(x), :)), 1), 2));
handles.ROI(handles.numTC).x  = x;
handles.ROI(handles.numTC).y = y;
handles.ROI(handles.numTC).z = z;
figure(2);
hold on;
plot(time, handles.ROI(handles.numTC).TC'+(10*handles.numTC), 'color', cmap(handles.numTC, :));

handles.numTC = handles.numTC+1;
guidata(hObject, handles);


% --- Executes on button press in clearTimecoursesPushbutton.
function clearTimecoursesPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to clearTimecoursesPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
choice = questdlg('Are you sure you want to delete all saved timecourses?', ...
	'Yes', 'No');
if (strcmp(choice, 'Yes'))
    try
        close(2)
    catch
    end
    handles.numTC = 1;
    handles.ROI = [];
end
guidata(hObject, handles);

% --- Executes on button press in exportTimecoursesPushbutton.
function exportTimecoursesPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to exportTimecoursesPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varname = cell2mat(inputdlg('Enter variable name of exported time courses:'));
assignin('base', varname, handles.ROI);

guidata(hObject, handles);
    


%%%
%%%  COLORMAP CALLBACKS
%%% 

function colorMaxSlider_Callback(hObject, eventdata, handles)
% hObject    handle to colorMaxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(handles.colorMaxSlider, 'Value');
set(handles.colorMaxTextbox, 'String', num2str(uint16(val*handles.maxDataColor)))
sliceslider_Callback(hObject, eventdata, handles)

function colorMinSlider_Callback(hObject, eventdata, handles)
% hObject    handle to colorMinSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(handles.colorMinSlider, 'Value');
set(handles.colorMinTextbox, 'String', num2str(uint16(val*handles.maxDataColor)))
sliceslider_Callback(hObject, eventdata, handles)

function colorMaxTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to colorMaxTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = str2num(get(handles.colorMaxTextbox, 'String'));
set(handles.colorMaxSlider, 'Value', val/handles.maxDataColor)
sliceslider_Callback(hObject, eventdata, handles)

function colorMinTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to colorMinTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = str2num(get(handles.colorMinTextbox, 'String'));
set(handles.colorMinSlider, 'Value', val/handles.maxDataColor)
sliceslider_Callback(hObject, eventdata, handles)

function logScaleCheckbox_Callback(hObject, eventdata, handles)
sliceslider_Callback(hObject, eventdata, handles)

%%%
%%% CROP FIELD CALLBACKS
%%%

function cropFieldCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to cropFieldCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(handles.cropFieldCheckbox, 'Value');
if(val == 0)
   set(handles.cropUIPanel, 'Visible', 'Off');
   handles.cropParams.xmin = handles.x(1);
   handles.cropParams.ymin = handles.y(1);
   handles.cropParams.zmin = handles.z(1);
   handles.cropParams.xmax = handles.x(end);
   handles.cropParams.ymax = handles.y(end);
   handles.cropParams.zmax = handles.z(end);
else
    set(handles.cropUIPanel, 'Visible', 'On');
    handles.cropParams.xmin = str2num(get(handles.cropXminTextbox, 'String'));
    handles.cropParams.ymin = str2num(get(handles.cropYminTextbox, 'String'));
    handles.cropParams.zmin = str2num(get(handles.cropZminTextbox, 'String'));
    handles.cropParams.xmax = str2num(get(handles.cropXmaxTextbox, 'String'));
    handles.cropParams.ymax = str2num(get(handles.cropYmaxTextbox, 'String'));
    handles.cropParams.zmax = str2num(get(handles.cropZmaxTextbox, 'String'));
end
sliceslider_Callback(hObject, eventdata, handles)

function cropXminTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to cropXminTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = str2num(get(handles.cropXminTextbox, 'String'));
if(val>=handles.cropParams.xmax)
    val = handles.cropParams.xmax-1;
    set(handles.cropXminTextbox, 'String', num2str(val));
elseif (val<0)
    val = 0;
    set(handles.cropXminTextbox, 'String', num2str(val));
end
cropFieldCheckbox_Callback(hObject, eventdata, handles)

function cropXmaxTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to cropXmaxTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = str2num(get(handles.cropXmaxTextbox, 'String'));
if(val<=handles.cropParams.xmin)
    val = handles.cropParams.xmin+1;
    set(handles.cropXmaxTextbox, 'String', num2str(val));
elseif (val>handles.x(end))
    val = handles.x(end);
    set(handles.cropXmaxTextbox, 'String', num2str(val));
end
cropFieldCheckbox_Callback(hObject, eventdata, handles)

function cropYminTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to cropYminTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = str2num(get(handles.cropYminTextbox, 'String'));
if(val>=handles.cropParams.ymax)
    val = handles.cropParams.ymax-1;
    set(handles.cropYminTextbox, 'String', num2str(val));
elseif (val<0)
    val = 0;
    set(handles.cropYminTextbox, 'String', num2str(val));
end
cropFieldCheckbox_Callback(hObject, eventdata, handles)

function cropYmaxTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to cropYmaxTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = str2num(get(handles.cropYmaxTextbox, 'String'));
if(val<=handles.cropParams.ymin)
    val = handles.cropParams.ymin+1;
    set(handles.cropYmaxTextbox, 'String', num2str(val));
elseif (val>handles.y(end))
    val = handles.y(end);
    set(handles.cropYmaxTextbox, 'String', num2str(val));
end
cropFieldCheckbox_Callback(hObject, eventdata, handles)

function cropZminTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to cropZminTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = str2num(get(handles.cropZminTextbox, 'String'));
if(val>=handles.cropParams.zmax)
    val = handles.cropParams.zmax-1;
    set(handles.cropZminTextbox, 'String', num2str(val));
elseif (val<0)
    val = 0;
    set(handles.cropZminTextbox, 'String', num2str(val));
end
cropFieldCheckbox_Callback(hObject, eventdata, handles)

function cropZmaxTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to cropZmaxTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = str2num(get(handles.cropZmaxTextbox, 'String'));
if(val<=handles.cropParams.zmin)
    val = handles.cropParams.zmin+1;
    set(handles.cropZmaxTextbox, 'String', num2str(val));
elseif (val>handles.z(end))
    val = handles.z(end);
    set(handles.cropZmaxTextbox, 'String', num2str(val));
end
cropFieldCheckbox_Callback(hObject, eventdata, handles)


%%%
%%% ASPECT RATIO CALLBACKS
%%%

% --- Executes on button press in uniformAspectRadiobutton.
function uniformAspectRadiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to uniformAspectRadiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of uniformAspectRadiobutton
set(handles.customCalPanel, 'Visible', 'Off');
handles.y = [1:handles.ss(1)]-1;
handles.x = [1:handles.ss(3)]-1;
handles.z = [1:handles.ss(2)]-1;
cropFieldCheckbox_Callback(hObject, eventdata, handles)

function calibratedAspectRadiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to calibratedAspectRadiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.y = ([1:handles.ss(1)]-1)*handles.info.GUIcalFactors.y_umPerPix;
try
    handles.x = ([1:handles.ss(3)]-1)*handles.info.GUIcalFactors.x_umPerPix;
catch
    handles.info.GUIcalFactors.x_umPerPix = handles.info.GUIcalFactors.xK_umPerVolt*handles.info.daq.scanAngle/(handles.info.daq.pixelsPerLine-2);
    handles.x = ([1:handles.ss(3)]-1)*handles.info.GUIcalFactors.x_umPerPix;
end
handles.z = ([1:handles.ss(2)]-1)*handles.info.GUIcalFactors.z_umPerPix;
set(handles.customCalPanel, 'Visible', 'Off');
cropFieldCheckbox_Callback(hObject, eventdata, handles)

function customAspectRadiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to customAspectRadiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.customCalPanel, 'Visible', 'On');

xCal = str2num(get(handles.xCalTextbox, 'String'));
yCal = str2num(get(handles.yCalTextbox, 'String'));
zCal = str2num(get(handles.zCalTextbox, 'String'));

handles.y = yCal*([1:handles.ss(1)]-1);
handles.x = xCal*([1:handles.ss(3)]-1);
handles.z = zCal*([1:handles.ss(2)]-1);
cropFieldCheckbox_Callback(hObject, eventdata, handles)

function xCalTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to xCalTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
customAspectRadiobutton_Callback(hObject, eventdata, handles)

function yCalTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to yCalTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
customAspectRadiobutton_Callback(hObject, eventdata, handles)

function zCalTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to zCalTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
customAspectRadiobutton_Callback(hObject, eventdata, handles)


%%%
%%% MIP CALLBACKS
%%%

function mipViewCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to mipViewCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(handles.mipViewCheckbox,'Value');
if (val)
    set(handles.mipUIPanel, 'Visible', 'On');
else
    set(handles.mipUIPanel, 'Visible', 'Off');
end
sliceslider_Callback(hObject, eventdata, handles)

function xSliceThicknessTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to xSliceThicknessTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% These 3 lines always make slice thickness an odd number (remember that
% going from -5 to 5 has 11 stops)

val = str2num(get(handles.xSliceThicknessTextbox, 'String'));
val = round((val-1)/2)*2+1;
val(val<0) = 0;
set(handles.xSliceThicknessTextbox, 'String', num2str(val));

sliceslider_Callback(hObject, eventdata, handles)

function ySliceThicknessTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to ySliceThicknessTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2num(get(handles.ySliceThicknessTextbox, 'String'));
val = round((val-1)/2)*2+1;
val(val<0) = 0;
set(handles.ySliceThicknessTextbox, 'String', num2str(val));

sliceslider_Callback(hObject, eventdata, handles)

function zSliceThicknessTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to zSliceThicknessTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2num(get(handles.zSliceThicknessTextbox, 'String'));
val = round((val-1)/2)*2+1;
val(val<0) = 0;
set(handles.zSliceThicknessTextbox, 'String', num2str(val));

sliceslider_Callback(hObject, eventdata, handles)


%%%
%%% CREATEFCNs GO AFTER THIS PORTION
%%%


% --- Executes during object creation, after setting all properties.
function colorMinTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorMinTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function timenum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timenum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function colorMaxSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorMaxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function slicenum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slicenum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function timeslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end





% --- Executes during object creation, after setting all properties.
function xyzview_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xyzview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function sliceslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes during object creation, after setting all properties.
function colorMinSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorMinSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes during object creation, after setting all properties.
function colorMaxTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorMaxTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function xCalTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xCalTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function yCalTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yCalTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function zCalTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zCalTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function cropXminTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cropXminTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function cropXmaxTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cropXmaxTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function cropYminTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cropYminTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function cropYmaxTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cropYmaxTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function cropZminTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cropZminTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function cropZmaxTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cropZmaxTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function xSliceThicknessTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xSliceThicknessTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function ySliceThicknessTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ySliceThicknessTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function zSliceThicknessTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zSliceThicknessTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.vid = [];
