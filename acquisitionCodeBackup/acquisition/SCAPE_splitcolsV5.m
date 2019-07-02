function varargout = SCAPE_splitcolsV5(varargin)
% SCAPE_SPLITCOLSV5 MATLAB code for SCAPE_splitcolsV5.fig
%      SCAPE_SPLITCOLSV5, by itself, creates a new SCAPE_SPLITCOLSV5 or
%      raises the existingparfor

%      singleton*.
%
%      H = SCAPE_SPLITCOLSV5 returns the handle to a new SCAPE_SPLITCOLSV5 or the handle to
%      the existing singleton*.
%
%      SCAPE_SPLITCOLSV5('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCAPE_SPLITCOLSV5.M with the given input arguments.
%
%      SCAPE_SPLITCOLSV5('Property','Value',...) creates a new SCAPE_SPLITCOLSV5 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SCAPE_splitcolsV5_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SCAPE_splitcolsV5_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SCAPE_splitcolsV5

% Last Modified by GUIDE v2.5 22-Aug-2016 16:35:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SCAPE_splitcolsV5_OpeningFcn, ...
    'gui_OutputFcn',  @SCAPE_splitcolsV5_OutputFcn, ...
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


% --- Executes just before SCAPE_splitcolsV5 is made visible.
function SCAPE_splitcolsV5_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SCAPE_splitcolsV5 (see VARARGIN)

% Choose default command line output for SCAPE_splitcolsV5
handles.output = hObject;
handles.inputdata = cell2mat(varargin(1));
handles.info = cell2mat(varargin(2));
ss = size(handles.inputdata);
infob{1} = sprintf('data size: %d (y) %d (z) %d (x)',ss(1), ss(2), ss(3));
if length(ss)==3
    ss(4) = 1;
end
infob{2} = sprintf('number of time-points: %d',ss(4));
handles.ss=ss;
handles.writetodesktop=0;

set(handles.infobox,'String',infob);

handles.viewstring = {'MIP x-y','MIP y-z'};
set(handles.MIPview,'String',handles.viewstring);
handles.justone=1;
figure(100)
clf
testtime = round(ss(4)/2);
temp = squeeze(max(handles.inputdata(:,:,:,testtime),[],3));
%  temp(1:floor(end/2),:) = 1*temp(1:floor(end/2),:);
imagesc(log(double(temp')));%imagesc(squeeze(max(moviesub_sh(1:end,:,:),[],1))');
title('select a region on each side that is the same, left to right')
colormap gray
[xg, zg] = ginput_col(1);
hold on
plot(xg,zg,'.r'); % EDIT Kripa - red and green switched for brain tissue
[xr, zr] = ginput_col(1);
hold on
plot(xr,zr,'.g');

clf
temp= squeeze(max(handles.inputdata(:,:,:,testtime),[],2));
%  temp(1:floor(end/2),:) = 5*temp(1:floor(end/2),:);
% temp(floor(end/2):end,:) = 5*temp(floor(end/2):end,:);
imagesc(log(double(temp')));%imagesc(squeeze(max(moviesub_sh(1:end,:,:),[],1))');
title('select a region on each side that is the same, left to right')
colormap gray
[xg, yg] = ginput_col(1);
hold on
plot(xg,yg,'.r'); % EDIT Kripa - red and green switched for brain tissue
[xr, yr] = ginput_col(1);
hold on
plot(xr,yr,'.g');
close
% handles.proc.rgb.xg2 = xg;
% handles.proc.rgb.xr2 = xr;
% handles.proc.rgb.yr2 = yr;
% handles.proc.rgb.yg2 = yg;
% handles.proc.rgb.zr2 = zr;
% handles.proc.rgb.zg2 = zg;

handles.proc.rgb.xg2 = 93.0660;
handles.proc.rgb.xr2 = 910.5090;
handles.proc.rgb.yr2 = 48.3025;
handles.proc.rgb.yg2 = 48.3025;
handles.proc.rgb.zr2 = 135.4537;
handles.proc.rgb.zg2 = 132.1636;

guidata(hObject,handles)
[red1, green1] = go_Callback(hObject, eventdata, handles,testtime);
set(handles.startframe,'String',num2str(testtime));
% check_Callback(hObject, eventdata, handles);

scr = 1/max(max(max(red1)));
scg = 1/max(max(max(green1)));
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

set(handles.scalered,'String',num2str(scr));%,'%.1f'));
set(handles.scalegreen,'String',num2str(scg));%,'%.1f'));
set(handles.minr,'String',num2str(minr));%,'%.1f'));
set(handles.ming,'String',num2str(ming));%,'%.1f'));

% subplot(1,2,1)
clear rgb
rgb(:,:,1) = squeeze(max(scr*(red1(:,:,3:end-3,1)-minr),[],2))';
rgb(:,:,2) = squeeze(scg*max((green1(:,:,3:end-3,1)-ming),[],2))';
rgb(:,:,3) = zeros(size(rgb(:,:,1)));
axes(handles.MainFigure)
imagesc(squeeze(uint8((256)*rgb)));
%  set(handles.viewstring,'Value',2)
set(handles.startframe,'String',num2str(testtime));
% Update handles structure
handles.writetodesktop=1;
handles.justone=0;

guidata(hObject, handles);

% UIWAIT makes SCAPE_splitcolsV5 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SCAPE_splitcolsV5_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function fileout_Callback(~, ~, ~)
% hObject    handle to fileout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileout as text
%        str2double(get(hObject,'String')) returns contents of fileout as a double


% --- Executes during object creation, after setting all properties.
function fileout_CreateFcn(hObject, ~, ~)
% hObject    handle to fileout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in go.
function [red, green] = go_Callback(hObject, ~, handles,timepoint)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Write to tiff
% profile on
% %%
xg = handles.proc.rgb.xg2;% = xg;
xr = handles.proc.rgb.xr2;% = xr;
yr = handles.proc.rgb.yr2;% = yr;
yg = handles.proc.rgb.yg2;% = yg;
zr = handles.proc.rgb.zr2;% = zr;
zg = handles.proc.rgb.zg2;% = zg;
%
handles.proc.rgb

if handles.writetodesktop ==1
    if length(size(handles.inputdata))==3
        timepoint= 1:1;
    else
        timepoint= 1:handles.ss(4);
    end
end

clear red green
left = floor(min([xg xr]))-1;
right = floor(size(handles.inputdata,1)-xr);

top = floor(min([zr zg]))-1;
bot = floor(size(handles.inputdata,2)-max([zr zg]));
topy = floor(min([yr yg]))-1;
boty = floor(size(handles.inputdata,3)-max([yr yg]));
% EDIT Kripa - red and green switched for brain tissue
if length(size(handles.inputdata))>3
%     keyboard
    red = double(handles.inputdata(round(xg+[-left:right]),round(zg+[-top:bot]),round(yg+[-topy:boty]),timepoint));
    green = double(handles.inputdata(round(xr+[-left:right]),round(zr+[-top:bot]),round(yr+[-topy:boty]),timepoint));
else
    red = double(handles.inputdata(round(xg+[-left:right]),round(zg+[-top:bot]),round(yg+[-topy:boty]),timepoint));
    green = double(handles.inputdata(round(xr+[-left:right]),round(zr+[-top:bot]),round(yr+[-topy:boty]),timepoint));
end

% rotation and tilt correction for imagesplitter
green = imrotate(green,str2num(get(handles.RotateGreen,'String')),'bilinear','Crop');
green = imresize(green,str2num(get(handles.ScaleGreen,'String')));

if str2num(get(handles.ScaleGreen,'String')) <1
    temp = zeros(size(red));
    temp(1:size(green,1),1:size(green,2),1:size(green,3))=green;
    green = temp;
else
    green = green(1:size(red,1),1:size(red,2),1:size(red,3),1:size(red,4));
end

if handles.writetodesktop ==1
    scr = str2num(get(handles.scalered,'String'));
    scg = str2num(get(handles.scalegreen,'String'));
    minr = str2num(get(handles.minr,'String'));
    ming = str2num(get(handles.ming,'String'));
    
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
    % temp1 = squeeze(handles.inputdata(:,:,:,i));
%     temp11 = squeeze(uint8(256/log(256)*log(256*rgb+1)));
%     temp12 = squeeze(uint8(256/log(256)*log(256*rgb2+1)));
%     figure(100);
%     clf
%     imagesc(temp11);
%     title('select crop area (lateral - top-left, then bottom-right)');
%     [y1, x1] = ginput_col(1);
%     hold on
%     plot(y1,x1,'.w');
%     [y2, x2] = ginput_col(1);
%     hold on
%     plot(y2,x2,'.w');
%     crp0 = floor(y1);
%     crp0y = floor(x1);
%     crp1 = ceil(y2);
%     crp1y = ceil(x2);
%     figure(100)
%     clf
%     imagesc(temp12)
%     title('select crop area (depth - right, then left)');
%     [z1, y1] = ginput_col(1);
%     hold on
%     plot(z1,y1,'.w');
%     hold on
%     [z2, y2] = ginput_col(1);
%     plot(z2,y2,'.w');
%     crp2 = floor(z1);
%     crp3 = ceil(z2);
%     close figure(100)
    
%     if length(size(red))==3
%         red = squeeze(double(red(crp0y:crp1y,crp2:crp3,crp0:crp1)));
%         green = squeeze(double(green(crp0y:crp1y,crp2:crp3,crp0:crp1)));
%     else
%         red = squeeze(double(red(crp0y:crp1y,crp2:crp3,crp0:crp1,timepoint)));
%         green = squeeze(double(green(crp0y:crp1y,crp2:crp3,crp0:crp1,timepoint)));
%     end
    
    assignin('base','red',red);
    assignin('base','green',green);
    assignin('base','params',handles.proc);
    assignin('base','scanName',handles.info.scanName);
end
handles.writetodesktop=1;
guidata(hObject,handles);
%
% pause
%
%

% profile viewer

% --- Executes on button press in check.
function check_Callback(hObject, eventdata, handles)
% hObject    handle to check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ss = size(handles.inputdata);
if length(handles.ss) == 3; handles.ss(4) = 1; end
startframe = str2num(get(handles.startframe,'String'));
testtime = startframe;
m=0;
if handles.justone==1; testtime = max(testtime); end
guidata(hObject,handles)
for i = 1:length(testtime)
    handles.writetodesktop=0;
    [red1, green1] = go_Callback(hObject, eventdata, handles, testtime(i));
    if get(handles.autoscale,'Value');
        scr = 1/max(max(max(red1)));
        scg = 1/max(max(max(green1)));
        minr = double(min(min(min(min(red1)))));
        ming = double(min(min(min(min(green1)))));
        scr =1/(1/scr-minr);
        scg =1/(1/scg-minr);
        handles.proc.autosc_scalered = scr;
        handles.proc.autosc_scalegreen = scg;
        handles.proc.autosc_minred = minr;
        handles.proc.autosc_mingreen = ming;
    else
        minr = str2num(get(handles.minr,'String'));%,'%.1f'));
        ming = str2num(get(handles.ming,'String'));%,'%.1f'));
        scr = str2num(get(handles.scalered,'String'));%,'%.1f'));
        scg = str2num(get(handles.scalegreen,'String'));%,'%.1f'));
     end
    handles.proc.scalered = scr;
    handles.proc.scalegreen = scg;
    handles.proc.minred = minr;
    handles.proc.mingreen = ming;
    set(handles.scalered,'String',num2str(scr));
    set(handles.scalegreen,'String',num2str(scg));
    set(handles.minr,'String',num2str(minr));
    set(handles.ming,'String',num2str(ming));
    
    if get(handles.MIPview,'Value')== 2;
        clear rgb
        rgb(:,:,1) = squeeze(max(scr*(red1(:,:,3:end-3,1)-minr),[],3)');
        rgb(:,:,2) = squeeze(scg*max((green1(:,:,3:end-3,1)-ming),[],3)');
        rgb(:,:,3) = zeros(size(rgb(:,:,1)));
        axes(handles.MainFigure)
        imagesc(squeeze(uint8((256)*rgb)));
    end
    if get(handles.MIPview,'Value')== 1;
        
        axis(handles.MainFigure);
        clear rgb
        rgb(:,:,1) = squeeze(max(scr*(red1(:,:,3:end-3,1)-minr),[],2))';
        rgb(:,:,2) = squeeze(scg*max(green1(:,:,3:end-3,1)-ming,[],2))';
        rgb(:,:,3) = zeros(size(rgb(:,:,1)));
        imagesc(squeeze(uint8((256)*rgb)));
        
    end
    title(['frame: ',num2str(testtime(i))])
    pause(0.1)
end
handles.justone = 0;
handles.writetodesktop=1;
guidata(hObject, handles)
%
%
% end

function startframe_Callback(~, ~, ~)
% hObject    handle to startframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startframe as text
%        str2double(get(hObject,'String')) returns contents of startframe as a double


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



% Hint: get(hObject,'Value') returns toggle state of yslices


% --- Executes on selection change in MIPview.
function MIPview_Callback(hObject, eventdata, handles)
% hObject    handle to MIPview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MIPview contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MIPview
check_Callback(hObject, eventdata, handles);


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

% Hints: get(hObject,'String') returns contents of interpval as text
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
% p = str2num(get(handles.startframe,'String'));
% set(handles.startframe,'String',num2str(p+1));
% if get(handles.MIPview,'Value')==1
handles.proc.rgb.xg2 = 1+ handles.proc.rgb.xg2;% = xg;
% end
% if get(handles.MIPview,'Value')==1
% handles.proc.rgb.xg2 = 1+ handles.proc.rgb.xg2;% = xg;
% end

% xr = handles.proc.rgb.xr2;% = xr;
% yr = handles.proc.rgb.yr2;% = yr;
% yg = handles.proc.rgb.yg2;% = yg;
% zr = handles.proc.rgb.zr2;% = zr;
% zg = handles.proc.rgb.zg2;% = zg;
% end
guidata(hObject,handles)
check_Callback(hObject, eventdata, handles);

% --- Executes on button press in minusleft.
function minusleft_Callback(hObject, eventdata, handles)
% hObject    handle to minusleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if get(handles.MIPview,'Value')==1
handles.proc.rgb.xg2 = -1+ handles.proc.rgb.xg2;% = xg;
% end
% if get(handles.MIPview,'Value')==1
% handles.proc.rgb.xg2 = 1+ handles.proc.rgb.xg2;% = xg;
% endp = str2num(get(handles.startframe,'String'));
% set(handles.startframe,'String',num2str(p-1));
guidata(hObject,handles)
check_Callback(hObject, eventdata, handles);


% --- Executes on button press in plustop.
function plustop_Callback(hObject, eventdata, handles)
% hObject    handle to plustop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.MIPview,'Value')==1
    handles.proc.rgb.yg2 = 1+ handles.proc.rgb.yg2;% = xg;
end
if get(handles.MIPview,'Value')==2
    handles.proc.rgb.zg2 = 1+ handles.proc.rgb.zg2;% = xg;
end
guidata(hObject,handles)
check_Callback(hObject, eventdata, handles);

% --- Executes on button press in minusbottom.
function minusbottom_Callback(hObject, eventdata, handles)
% hObject    handle to minusbottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.MIPview,'Value')==1
    handles.proc.rgb.yg2 = -1+ handles.proc.rgb.yg2;% = xg;
end
if get(handles.MIPview,'Value')==2
    handles.proc.rgb.zg2 = -1+ handles.proc.rgb.zg2;% = xg;
end
guidata(hObject,handles)
check_Callback(hObject, eventdata, handles);





function scalered_Callback(hObject, ~, handles)
% hObject    handle to scalered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scalered as text
%        str2double(get(hObject,'String')) returns contents of scalered as a double
handles.proc.scalered = str2num(get(handles.scalered,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function scalered_CreateFcn(hObject, ~, ~)
% hObject    handle to scalered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scalegreen_Callback(hObject, eventdata, handles)
% hObject    handle to scalegreen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scalegreen as text
%        str2double(get(hObject,'String')) returns contents of scalegreen as a double
handles.proc.scalegreen = str2num(get(handles.scalegreen,'String'));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function scalegreen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scalegreen (see GCBO)
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
% if isfield(handles,'info');
% [outdir] = handles.info.dataDirectory;
% else
[outdir] = uigetdir('Select output directory');
% end
dataDirectory = outdir;
if(isdir([dataDirectory , '/', 'tiff_stacks']) == 0)
    mkdir ([dataDirectory , '/','tiff_stacks']);
end
scanName = handles.info.scanName;

if (isdir([dataDirectory '/tiff_stacks/', scanName]) == 0)
    mkdir([dataDirectory '/tiff_stacks/', scanName]);
end
%             vals = str2num(get(handles.interpval,'String'));

filepath = [dataDirectory '/tiff_stacks/', scanName,'/'];

% need to decide whether to use the ones for the current frame (auto or user), or over the
% whole dataset

scr = str2num(get(handles.scalered,'String'));
scg = str2num(get(handles.scalegreen,'String'));
minr = str2num(get(handles.minr,'String'));
ming = str2num(get(handles.ming,'String'));

handles.writetodesktop=0;
guidata(hObject,handles);
if get(handles.preview_crop,'Value')
    % quick check
    i = max(str2num(get(handles.startframe,'String')));
    [red, green] = go_Callback(hObject, eventdata, handles,i);
    chk = 'replay';
    while strcmp(chk,'replay');
        figure(100);
        clf
        clear rgb
        for j = fliplr([round(size(red,2)/10):round(size(red,2)/10):size(red,2)-round(size(red,2)/10)]);
            rgb(:,:,1) = squeeze(scr*(double(red(:,j,:))-minr));
            rgb(:,:,2) = squeeze(scg*(double(green(:,j,:))-ming));
            rgb(:,:,3) = zeros(size(rgb(:,:,1)));
            % temp1 = squeeze(handles.inputdata(:,:,:,i));
            temp11 = squeeze(uint8(256*rgb));%(2^16)*temp1/(max(max(max(temp1))));
            
            for i = 1:2
                subplot(1,3,i)
                imagesc(squeeze(temp11(:,:,i)))
                colorbar
                colormap jet
                axis tight
                
            end
            subplot(1,3,3)
            imagesc(temp11);
            axis tight
            pause(0.1)
            
        end
        chk = questdlg('look ok?','prepare for tifs','replay','crop','yes','no');
    end
    chk = 'crop';
    if strcmp(chk,'crop')
        close figure 100
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
        % temp1 = squeeze(handles.inputdata(:,:,:,i));
        temp11 = squeeze(uint8(256*rgb));
        temp12 = squeeze(uint8(256*rgb2));
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
        title('select crop area (depth - right, then left)');
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
    else
        crp0 = 1;
        crp1 = size(red,3);
        crp0y = 1;
        crp1y = size(red,1);
        crp2 = 1;
        crp3 = size(red,2);
    end
else
    chk='yes';
    crp0 = 0;
    crp1 = 0;
    crp0y = 0;
    crp1y = 0;
    crp2 = 0;
    crp3 = 0;
end
if strcmp(chk,'yes')
    disp('Writing tiffs. Please wait.')
    m=0;
    conversionFactors = [handles.info.cal.ylat; handles.info.cal.zdep; handles.info.cal.xwid];
    if length(handles.ss)==3; handles.ss(4) = 1; end
    if handles.ss(4)==1
        if get(handles.SeparateTiffs,'Value')
            imgToSave1 = [filepath,'RGB_' scanName, '_R.tiff'];
            imgToSave2 = [filepath,'RGB_' scanName, '_G.tiff'];
        else
            imgToSave = [filepath,'uncorrected_RGB_' scanName, '.tiff'];
        end
        %normalize each volume to its max
        [red, green] = go_Callback(hObject, eventdata, handles,1);
        if crp0==0; crp0 = 1; crp1 = size(red,3); crp0y = 1; crp1y = size(red,1); crp2 = 1; crp3 = size(red,2); end
        if get(handles.SeparateTiffs,'Value')
            red = squeeze((double(red(crp0y:crp1y,crp2:crp3,crp0:crp1))));
            green = squeeze((double(green(crp0y:crp1y,crp2:crp3,crp0:crp1))));
        else
            red = squeeze(scr*(double(red(crp0y:crp1y,crp2:crp3,crp0:crp1))-minr));
            green = squeeze(scg*(double(green(crp0y:crp1y,crp2:crp3,crp0:crp1))-ming));
        end
        if get(handles.Skewcorrection,'Value')
            if get(handles.SeparateTiffs,'Value')
                imgToSave1 = [filepath,'RGB_' scanName, '_R.tiff'];
                imgToSave2 = [filepath,'RGB_' scanName, '_G.tiff'];
            else
                imgToSave = [filepath,'RGB_' scanName, '.tiff'];
            end
            % Coordinate System Correction
            red = flip(flip(red,2),3);
            red = permute(red, [3 1 2]);
            RR = imref3d(size(red), conversionFactors(3), conversionFactors(1), conversionFactors(2));
            
            green = flip(flip(green(:,:,:),2),3);
            green = permute(green, [3 1 2]);
            
            RG = imref3d(size(green), conversionFactors(3), conversionFactors(1), conversionFactors(2));
            % Correct for skew
%             
%             figure(99); subplot(2,1,1)
%             imagesc([0:size(red, 1)-1]*conversionFactors(3),[0:size(red, 3)-1]*conversionFactors(2), log(squeeze(max(red(:,round(end/3):round(2*end/3), :), [], 2)))')
%             colormap gray; axis image; title('Raw')
            
            delta = str2num(get(handles.SkewAng,'String'));
            %delta = pi/2;
            affineMatrix = [1 0 0 0;
                0 1 0 0;
                0 cotd(delta) 1 0;
                0 0 0 1];
            tform = affine3d(affineMatrix);
            
            [red, ~] = imwarp(red, RR, tform);
            [green, ~] = imwarp(green, RG, tform);
%            figure(99); subplot(2,1,2)
% imagesc([0:size(red, 1)-1]*conversionFactors(3),[0:size(red, 3)-1]*conversionFactors(2), log(squeeze(max(red(:,round(end/3):round(2*end/3), :), [], 2)))')
% colormap gray; axis image;  title('Skew Corrected')
%             disp('Skew correction applied')
        end
        ss = size(red);
        for i = 1:ss(3)
            if get(handles.SeparateTiffs,'Value')
                imwrite(uint16(red(:,:,i)), imgToSave1, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                imwrite(uint16(green(:,:,i)), imgToSave2, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
            else
                clear rgb
                rgb(:,:,1) = red(:,:,i);
                rgb(:,:,2) = green(:,:,i);
                rgb(:,:,3) = zeros(size(red(:,:,i)));
                temp11 = squeeze(uint16(2^16*rgb));%(2^16)*temp1/(max(max(max(temp1))));
                
                imwrite(temp11, imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
            end
        end
        
        if get(handles.MIP, 'Value')
            filename =  fullfile(dataDirectory, ['RGB' scanName '__HIRES.tif']);
            tempr = squeeze(max(red(:,:,:),[],2))';
            tempr = cat(1,tempr,1000*ones(2,size(red,1)));
            tempr = cat(1,tempr, (squeeze(max(red,[],3))'));
            tempg = squeeze(max(green(:,:,:),[],2))';
            tempg = cat(1,tempg,1000*ones(2,size(green,1)));
            tempg = cat(1,tempg, (squeeze(max(green,[],3))'));
            temprg(:,:,1) = tempr;
            temprg(:,:,2) = tempg;
            temprg(:,:,3) = zeros(size(tempr));
            imwrite(temprg,filename,'tif', 'Compression', 'none')
        end
        
    else
        [red, ~] = go_Callback(hObject, eventdata, handles,1);
        if crp0==0; crp0 = 1; crp1 = size(red,3); crp0y = 1; crp1y = size(red,1); crp2 = 1; crp3 = size(red,2); end
        
%         if get(handles.MIP,'Value')
%             filename =  fullfile(dataDirectory, ['RGB' scanName '__movie.avi']);
%             vidObj = VideoWriter(filename, 'Uncompressed AVI');
%             vidObj.FrameRate = handles.info.daq.scanRate/handles.info.daq.scanLength;
%             open(vidObj);
%         end
        
        
        parfor_progress(handles.ss(4));
        tic
        for i = 1:handles.ss(4) %scans
            if get(handles.SeparateTiffs,'Value')
                imgToSave1 = [filepath,'uncorrected_RGB_' scanName, '_', num2str(i) '_R.tiff'];
                imgToSave2 = [filepath,'uncorrected_RGB_' scanName, '_', num2str(i) '_G.tiff'];
            else
                imgToSave = [filepath,'uncorrected_RGB_' scanName, '_', num2str(i) '.tiff'];
            end
            %normalize each volume to its max
            [red, green] = go_Callback(hObject, eventdata, handles,i);
            red = squeeze((double(red(crp0y:crp1y,crp2:crp3,crp0:crp1))-0));
            green = squeeze((double(green(crp0y:crp1y,crp2:crp3,crp0:crp1))-0));
            
            
            if get(handles.Skewcorrection,'Value')
                if get(handles.SeparateTiffs,'Value')
                    imgToSave1 = [filepath,'G_' scanName, '_', num2str(i) '.tiff'];
                    imgToSave2 = [filepath,'R_' scanName, '_', num2str(i) '.tiff'];
                else
                    imgToSave = [filepath,'RGB_' scanName, '_', num2str(i) '.tiff'];
                end
                % Coordinate System Correction
                red = flip(flip(red(:,:,:),2),3);
                red = permute(red, [3 1 2]);
                RR = imref3d(size(red), conversionFactors(3), conversionFactors(1), conversionFactors(2));
                
                green = flip(flip(green(:,:,:),2),3);
                green = permute(green, [3 1 2]);
                RG = imref3d(size(green), conversionFactors(3), conversionFactors(1), conversionFactors(2));
                
                % Correct for skew
                %                     figure(99); subplot(1,2,1)
                %                     imagesc([0:size(red, 1)-1]*conversionFactors(3),[0:size(red, 3)-1]*conversionFactors(2), squeeze(max(red(:,round(end/3):round(2*end/3), :), [], 2))')
                %                     colormap gray; axis image
                
                delta = str2num(get(handles.SkewAng,'String'));
                %delta = pi/2;
                affineMatrix = [1 0 0 0;
                    0 1 0 0;
                    0 cotd(delta) 1 0;
                    0 0 0 1];
                tform = affine3d(affineMatrix);
                [red , ~] = imwarp(red, RR, tform);
                [green , ~] = imwarp(green, RG, tform);
            end
            ss = size(red);
            
            for j = 1:ss(3)
                if get(handles.SeparateTiffs,'Value')
                    imwrite(uint16(red(:,:,j)), imgToSave1, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                    imwrite(uint16(green(:,:,j)), imgToSave2, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                else
                    rgb = cat(3, red(:, :, j), green(:, :, j), zeros(size(red(:,:,j))));
                    imwrite(squeeze(uint16(2^16*rgb)), imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                    
                end
            end
            
%             if get(handles.MIP,'Value')
%                 tmpr = max(max(red(:,:,:),[],3)-min(red(:,:,:),[],3),[],2);
%                 [ar, br] = hist(reshape(double(tmpr),[1,prod(size(tmpr))]),200);
%                 cc1r = br(min(find(ar>0.01*mean(ar))));
%                 cc2r = br(max(find(ar>0.01*mean(ar))));
%                 tmpg = max(max(green(:,:,:),[],3)-min(green(:,:,:),[],3),[],2);
%                 [ag, bg] = hist(reshape(double(tmpg),[1,prod(size(tmpg))]),200);
%                 cc1g = bg(min(find(ag>0.01*mean(ag))));
%                 cc2g = bg(max(find(ag>0.01*mean(ag))));
%                 
%                 tempr = squeeze(max(red(:,2:end-2,:)-red(:,:,:),[],2))';
%                 tempr = cat(1,tempr,1000*ones(2,size(red,1)));
%                 tempr = cat(1,tempr, (squeeze(max(red(:,2:end-2,:)-red(:,:,:),[],3))'));
%                 tempr = uint8((256/(cc2r-cc1r))*(tempr-cc1r));
%                 tempg = squeeze(max(green(:,2:end-2,:)-green(:,:,:),[],2))';
%                 tempg = cat(1,tempg,1000*ones(2,size(green,1)));
%                 tempg = cat(1,tempg, (squeeze(max(green(:,2:end-2,:)-green(:,:,:),[],3))'));
%                 tempg = uint8((256/(cc2g-cc1g))*(tempg-cc1g));
%                 temprg = cat(3,tempr, tempg, zeros(size(tempr)));
%                 M(i) = im2frame(temprg);
%             end
            
           parfor_progress;
        end
        if get(handles.MIP, 'Value')
            writeVideo(vidObj,M);
            close(vidObj)
        end
        if get(handles.Skewcorrection,'Value')
            %                     figure(99); subplot(1,2,2)
            %                     imagesc([0:size(red, 1)-1]*conversionFactors(3),[0:size(red, 3)-1]*conversionFactors(2), squeeze(max(red(:,round(end/3):round(2*end/3), :,1), [], 2))')
            %                     colormap gray; axis image
            
            disp('Skew correction applied')
        end
       parfor_progress(0);
        toc
    end
    disp('Finished writing tiffs');
else
    disp('cancelling')
end



% --- Executes on slider movement.
function redslide_Callback(hObject, eventdata, handles)
% hObject    handle to redslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
scr = handles.proc.scalered;%
s = 1- get(handles.redslide,'Value');
range = [0:2/scr];
scr = 1/(s*2/scr);
set(handles.scalered,'String',num2str(scr));
set(handles.autoscale,'Value',0);
handles.justone=1;
guidata(hObject,handles)
check_Callback(hObject, eventdata, handles);
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
s = 1- get(handles.greenslide,'Value');
scg = 1/(s*2/scg);
set(handles.scalegreen,'String',num2str(scg));
set(handles.autoscale,'Value',0);
handles.justone=1;
guidata(hObject,handles)
check_Callback(hObject, eventdata, handles);
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

if get(handles.autoscale,'Value');
    set(handles.scalegreen,'String',handles.proc.autosc_scalegreen);
    set(handles.greenslide,'Value',0.5);
    handles.proc.scalegreen = handles.proc.autosc_scalegreen;
    set(handles.scalered,'String',handles.proc.autosc_scalered);
    set(handles.redslide,'Value',0.5);
    handles.proc.scalered = handles.proc.autosc_scalered;
    set(handles.minred,'String',handles.proc.autosc_minred);
    set(handles.mingreen,'String',handles.proc.autosc_minred);
    check_Callback(hObject, eventdata, handles);
    guidata(hObject,handles)
else
end
% Hint: get(hObject,'Value') returns toggle state of autoscale



function minr_Callback(hObject, eventdata, handles)
% hObject    handle to minr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minr as text
%        str2double(get(hObject,'String')) returns contents of minr as a double
handles.justone=1;
guidata(hObject,handles)
check_Callback(hObject, eventdata, handles);
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

% Hints: get(hObject,'String') returns contents of ming as text
%        str2double(get(hObject,'String')) returns contents of ming as a double
handles.justone=1;
guidata(hObject,handles)
check_Callback(hObject, eventdata, handles);
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

% Hints: get(hObject,'String') returns contents of RotateGreen as text
%        str2double(get(hObject,'String')) returns contents of RotateGreen as a double


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



function ScaleGreen_Callback(hObject, eventdata, handles)
% hObject    handle to ScaleGreen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScaleGreen as text
%        str2double(get(hObject,'String')) returns contents of ScaleGreen as a double


% --- Executes during object creation, after setting all properties.
function ScaleGreen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScaleGreen (see GCBO)
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
[red, ~] = go_Callback(hObject, eventdata, handles,floor(handles.ss(4)/2));
% Coordinate System Correction
conversionFactors = [handles.info.cal.ylat; handles.info.cal.zdep; handles.info.cal.xwid];

red = flip(flip(red,2),3);
red = permute(red, [3 1 2]);
RR = imref3d(size(red), conversionFactors(3), conversionFactors(1), conversionFactors(2));

% Correct for skew

figure(99); subplot(2,1,1)
imagesc([0:size(red, 1)-1]*conversionFactors(3),[0:size(red, 3)-1]*conversionFactors(2), log(squeeze(max(red(:,round(end/3):round(2*end/3), :), [], 2)))')
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
imagesc([0:size(red, 1)-1]*conversionFactors(3),[0:size(red, 3)-1]*conversionFactors(2), log(squeeze(max(red(:,round(end/3):round(2*end/3), :), [], 2)))')
colormap gray; axis image;  title('Skew Corrected')
disp('Skew correction applied')
