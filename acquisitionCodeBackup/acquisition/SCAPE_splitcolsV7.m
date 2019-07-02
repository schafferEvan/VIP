function varargout = SCAPE_splitcolsV7(varargin)
% SCAPE_SPLITCOLSV7 MATLAB code for SCAPE_splitcolsV7.fig
%      SCAPE_SPLITCOLSV7, by itself, creates a new SCAPE_SPLITCOLSV7 or
%      raises the existingparfor
%      singleton*.
%
%      H = SCAPE_SPLITCOLSV7 returns the handle to a new SCAPE_SPLITCOLSV7 or the handle to
%      the existing singleton*.
%
%      SCAPE_SPLITCOLSV7('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCAPE_SPLITCOLSV7.M with the given input arguments.
%
%      SCAPE_SPLITCOLSV7('Property','Value',...) creates a new SCAPE_SPLITCOLSV7 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SCAPE_splitcolsV7_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SCAPE_splitcolsV7_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SCAPE_splitcolsV7

% Last Modified by GUIDE v2.5 31-Jul-2017 09:55:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SCAPE_splitcolsV7_OpeningFcn, ...
    'gui_OutputFcn',  @SCAPE_splitcolsV7_OutputFcn, ...
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


% --- Executes just before SCAPE_splitcolsV7 is made visible.
function SCAPE_splitcolsV7_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SCAPE_splitcolsV7 (see VARARGIN)

% Choose default command line output for SCAPE_splitcolsV7
handles.output = hObject;
handles.inputdata = cell2mat(varargin(1));
handles.info = cell2mat(varargin(2));
handles.dataDirectory = cell2mat(varargin(3));
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


% initial values for rotation and scaling
handles.scalepastval = str2num(get(handles.ScaleG,'String'));
handles.rotpastval = str2num(get(handles.RotateGreen,'String'));


if(isdir([handles.dataDirectory , '/', 'tiff_stacks']) == 0)
    mkdir ([handles.dataDirectory , '/','tiff_stacks']);
end
FilePath = dir(fullfile(handles.dataDirectory, 'tiff_stacks/*.mat'));
FilePath.name
if numel(FilePath)~=0
    presentDirectory = pwd;
    cd(fullfile(handles.dataDirectory, 'tiff_stacks/'));
    [fileName, pathName, ~] = uigetfile('*.mat');
    cd (pathName)
    transforms = load(fileName);
    cd(presentDirectory)
    
    set(handles.scalered,'String',num2str(transforms.transforms.scr));
    set(handles.scalegreen2,'String',num2str(transforms.transforms.scg));
    set(handles.minr,'String',num2str(transforms.transforms.minr));
    set(handles.ming,'String',num2str(transforms.transforms.ming));
    handles.proc.scalered = transforms.transforms.scr;
    handles.proc.scalegreen = transforms.transforms.scg;
    handles.proc.autosc_scalered = transforms.transforms.scr;
    handles.proc.autosc_scalegreen = transforms.transforms.scg;
    handles.proc.autosc_minred = transforms.transforms.minr;
    handles.proc.autosc_mingreen = transforms.transforms.ming;
    handles.proc.autosc_minred = transforms.transforms.minr;
    handles.proc.autosc_mingreen = transforms.transforms.ming;
    
    handles.proc.rgb.xg2 = transforms.transforms.xg;
    handles.proc.rgb.xr2 = transforms.transforms.xr;
    handles.proc.rgb.yr2 = transforms.transforms.yr;
    handles.proc.rgb.yg2 = transforms.transforms.yg;
    handles.proc.rgb.zr2 = transforms.transforms.zr;
    handles.proc.rgb.zg2 = transforms.transforms.zg;
    set(handles.RotateGreen,'String',num2str(transforms.transforms.rotateg));
    set(handles.ScaleG,'String',num2str(transforms.transforms.scaleg));
    
    testtime = 1;
    
    guidata(hObject,handles)
    tic
    clear red1 green1
    [red1, green1] = go_Callback(hObject, eventdata, handles, testtime);
    toc
    set(handles.startframe,'String',num2str(testtime));
    check_Callback(hObject, eventdata, handles);
    
    
    scr = str2num(get(handles.scalered,'String'));
    scg = str2num(get(handles.scalegreen2,'String'));
    minr = str2num(get(handles.minr,'String'));
    ming = str2num(get(handles.ming,'String'));
    
    clear rgb
    rgb(:,:,1) = squeeze(max(scr*double(red1(:,:,3:end-3,1)-minr),[],2))';
    rgb(:,:,2) = squeeze(scg*(max(double(green1(:,:,3:end-3,1))-ming,[],2)))';
    rgb(:,:,3) = zeros(size(rgb(:,:,1)));
    axis(handles.MainFigure)
    if get(handles.takelog,'Value')
        rgb = log(double(rgb).*256+1)./log(256);
    end
    imagesc(squeeze(uint8((256)*rgb)));
    handles.RGBMIPS = rgb;
    axis off
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
    plot(xg,zg,'.g'); % EDIT Kripa - red and green switched for homebuilt imagesplitter
    [xr, zr] = ginput_col(1);
    hold on
    plot(xr,zr,'.r');
    
    clf
    temp= squeeze(max(handles.inputdata(:,:,:,testtime),[],2));
    %  temp(1:floor(end/2),:) = 5*temp(1:floor(end/2),:);
    % temp(floor(end/2):end,:) = 5*temp(floor(end/2):end,:);
    imagesc(log(double(temp))');%imagesc(squeeze(max(moviesub_sh(1:end,:,:),[],1))');
    title('select a region on each side that is the same, left to right')
    colormap gray
    [xg, yg] = ginput_col(1);
    hold on
    plot(xg,yg,'.g'); % EDIT Kripa - red and green switched for brain tissue
    [xr, yr] = ginput_col(1);
    hold on
    plot(xr,yr,'.r');
    close
    handles.proc.rgb.xg2 = xg;
    handles.proc.rgb.xr2 = xr;
    handles.proc.rgb.yr2 = yr;
    handles.proc.rgb.yg2 = yg;
    handles.proc.rgb.zr2 = zr;
    handles.proc.rgb.zg2 = zg;
    guidata(hObject,handles)
    clear red1 green1
    [green1, red1] = go_Callback(hObject, eventdata, handles,testtime);
    set(handles.startframe,'String',num2str(testtime));
    check_Callback(hObject, eventdata, handles);
    
    
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
    set(handles.scalegreen2,'String',num2str(scg));
    set(handles.minr,'String',num2str(minr));
    set(handles.ming,'String',num2str(ming));
    
    % subplot(1,2,1)
    clear rgb
    rgb(:,:,1) = squeeze(max(scr*(double(red1(:,:,3:end-3,1))-minr),[],2))';
    rgb(:,:,2) = squeeze(scg*max((double(green1(:,:,3:end-3,1))-ming),[],2))';
    rgb(:,:,3) = zeros(size(rgb(:,:,1)));
    if get(handles.takelog,'Value')
        rgb = log(double(rgb).*256+1)./log(256);
    end
    axis(handles.MainFigure)
    imagesc(squeeze(uint8((256)*rgb)));
    axis off
    handles.RGBMIPS = rgb;
    
    %  set(handles.viewstring,'Value',2)
    set(handles.startframe,'String',num2str(testtime));
    % Update handles structure
    handles.writetodesktop=1;
    handles.justone=0;
end

guidata(hObject, handles);

% UIWAIT makes SCAPE_splitcolsV7 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SCAPE_splitcolsV7_OutputFcn(~, ~, handles)
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

xg = handles.proc.rgb.xg2;
xr = handles.proc.rgb.xr2;
yr = handles.proc.rgb.yr2;
yg = handles.proc.rgb.yg2;
zr = handles.proc.rgb.zr2;
zg = handles.proc.rgb.zg2;

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
% EDIT Kripa - red and green set for new imagesplitter
if length(size(handles.inputdata))>3
    green = (handles.inputdata(round(xg+[-left:right]),round(zg+[-top:bot]),round(yg+[-topy:boty]),timepoint));
    red = (handles.inputdata(round(xr+[-left:right]),round(zr+[-top:bot]),round(yr+[-topy:boty]),timepoint));
else
    green = (handles.inputdata(round(xg+[-left:right]),round(zg+[-top:bot]),round(yg+[-topy:boty])));
    red = (handles.inputdata(round(xr+[-left:right]),round(zr+[-top:bot]),round(yr+[-topy:boty])));
end

% rotation correction for 1st gen imagesplitter
if (str2num(get(handles.RotateGreen,'String'))~=0)
    green = imrotate(green,str2num(get(handles.RotateGreen,'String')),'bilinear','Crop');
end
% scale correction for 1st gen imagesplitter


if (str2num(get(handles.ScaleG,'String'))~=1)
    green = imresize(green,str2num(get(handles.ScaleG,'String')));
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
    clear red1 green1
    [red1, green1] = go_Callback(hObject, eventdata, handles, testtime(i));
    %     if get(handles.autoscale,'Value');
    
    minr = str2num(get(handles.minr,'String'));%,'%.1f'));
    ming = str2num(get(handles.ming,'String'));%,'%.1f'));
    scr = str2num(get(handles.scalered,'String'));%,'%.1f'));
    scg = str2num(get(handles.scalegreen2,'String'));%,'%.1f'));
    
    handles.proc.scalered = scr;
    handles.proc.scalegreen = scg;
    handles.proc.minred = minr;
    handles.proc.mingreen = ming;
    
    if get(handles.MIPview,'Value')== 2;
        clear rgb
        rgb(:,:,1) = squeeze(max(scr*(double(red1(:,:,round(end/5):round(4*end/5),1))-minr),[],3)');
        rgb(:,:,2) = squeeze(scg*max((double(green1(:,:,round(end/5):round(4*end/5),1))-ming),[],3)');
        rgb(:,:,3) = zeros(size(rgb(:,:,1)));
        if get(handles.takelog,'Value')
            rgb = log(double(rgb).*256+1)./log(256);
        end
        axis(handles.MainFigure)
        imagesc(squeeze(uint8((256)*rgb)));
        handles.RGBMIPS = rgb;
    end
    if get(handles.MIPview,'Value')== 1;
        
        axis(handles.MainFigure);
        clear rgb
        rgb(:,:,1) = squeeze(max(scr*(double(red1(:,round(end/5):round(4*end/5),3:end-3,1))-minr),[],2))';
        rgb(:,:,2) = squeeze(scg*max(double(green1(:,round(end/5):round(4*end/5),3:end-3,1))-ming,[],2))';
        rgb(:,:,3) = zeros(size(rgb(:,:,1)));
        if get(handles.takelog,'Value')
            rgb = log(double(rgb).*256+1)./log(256);
        end
        max(max(max(rgb)))
        min(min(min(rgb)))
        size(rgb)
        imagesc(uint8((256)*rgb));
        handles.RGBMIPS = rgb;
    end
    title(['frame: ',num2str(testtime(i))])
    pause(0.1)
end
handles.justone = 0;
handles.writetodesktop=1;
guidata(hObject, handles)


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

% if get(handles.MIPview,'Value')==1
handles.proc.rgb.xg2 = -1+ handles.proc.rgb.xg2; % = xg;
% end
% if get(handles.MIPview,'Value')==1
% handles.proc.rgb.xg2 = 1+ handles.proc.rgb.xg2;% = xg;
% end

axis(handles.MainFigure);
clear rgb
rgb(:,:,1) = padarray(handles.RGBMIPS(2:end-1,2:end-1,1), [1,1]);
rgb(:,:,2) = padarray(handles.RGBMIPS(2:end-1,1:end-2,2), [1,1]);
rgb(:,:,3) = zeros(size(rgb(:,:,1)));
imagesc(squeeze(uint8((256)*rgb)));
handles.RGBMIPS = rgb;
guidata(hObject,handles)
display('image adjusted')
%check_Callback(hObject, eventdata, handles);

% --- Executes on button press in minusleft.
function minusleft_Callback(hObject, eventdata, handles)
% hObject    handle to minusleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if get(handles.MIPview,'Value')==1
handles.proc.rgb.xg2 = 1+ handles.proc.rgb.xg2; % = xg;
% end
% if get(handles.MIPview,'Value')==1
% handles.proc.rgb.xg2 = 1+ handles.proc.rgb.xg2;% = xg;
% endp = str2num(get(handles.startframe,'String'));
% set(handles.startframe,'String',num2str(p-1));

%check_Callback(hObject, eventdata, handles);
axis(handles.MainFigure);
clear rgb
rgb(:,:,1) = padarray(handles.RGBMIPS(2:end-1,2:end-1,1),[1,1]);
rgb(:,:,2) = padarray(handles.RGBMIPS(2:end-1,3:end,2),[1,1]);
rgb(:,:,3) = zeros(size(rgb(:,:,1)));
imagesc(squeeze(uint8((256)*rgb)));
handles.RGBMIPS = rgb;
guidata(hObject,handles)
display('image adjusted')


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
axis(handles.MainFigure);
clear rgb
rgb(:,:,1) = padarray(handles.RGBMIPS(2:end-1,2:end-1,1),[1,1]);
rgb(:,:,2) = padarray(handles.RGBMIPS(3:end,2:end-1,2),[1,1]);
rgb(:,:,3) = zeros(size(rgb(:,:,1)));
imagesc(squeeze(uint8((256)*rgb)));
handles.RGBMIPS = rgb;
guidata(hObject,handles)
display('image adjusted')


%check_Callback(hObject, eventdata, handles);

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
axis(handles.MainFigure);
clear rgb
rgb(:,:,1) = padarray(handles.RGBMIPS(2:end-1,2:end-1,1),[1,1]);
rgb(:,:,2) = padarray(handles.RGBMIPS(1:end-2,2:end-1,2),[1,1]);
rgb(:,:,3) = zeros(size(rgb(:,:,1)));
imagesc(squeeze(uint8((256)*rgb)));

handles.RGBMIPS=rgb;
guidata(hObject,handles)
display('image adjusted')

%check_Callback(hObject, eventdata, handles);


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




% --- Executes on button press in tiffs.
function tiffs_Callback(hObject, eventdata, handles)
% hObject    handle to tiffs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% scans = str2num(get(handles.startframe,'String'));
% ss = size(handles.inputdata);
% % need to handle normalization better
if isfield(handles,'info');
    [outdir] = handles.dataDirectory;
else
    [outdir] = uigetdir('Select output directory');
end

dataDirectory = outdir;
if(isdir([dataDirectory , '/', 'tiff_stacks']) == 0)
    mkdir ([dataDirectory , '/','tiff_stacks']);
end
scanName = handles.info.scanName;

if (isdir([dataDirectory '/tiff_stacks/', scanName]) == 0)
    mkdir([dataDirectory '/tiff_stacks/', scanName]);
end
%             vals = str2num(get(handles.interpval,'String'));

FilePath = [dataDirectory '/tiff_stacks/', scanName,'/'];

% need to decide whether to use the ones for the current frame (auto or user), or over the
% whole dataset

scr = str2num(get(handles.scalered,'String'));
scg = str2num(get(handles.scalegreen2,'String'));
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
            rgb(:,:,1) = squeeze(scr.*(double(red(:,j,:))-minr));
            rgb(:,:,2) = squeeze(scg.*(double(green(:,j,:))-ming));
            rgb(:,:,3) = zeros(size(rgb(:,:,1)));
            if get(handles.takelog,'Value')
                rgb = log(double(rgb).*256+1)./log(256);
            end
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
            if get(handles.takelog,'Value')
                rgb = log(double(rgb).*256+1)./log(256);
            end
            rgb2(:,:,1) = squeeze(scr*(double(max(red(:,:,:),[],3))-minr));
            rgb2(:,:,2) = squeeze(scg*(double(max(green(:,:,:),[],3))-ming));
            rgb2(:,:,3) = zeros(size(rgb2(:,:,1)));
            if get(handles.takelog,'Value')
                rgb2 = log(double(rgb2).*256+1)./log(256);
            end
        else
            rgb(:,:,1) = squeeze(scr*(double(max(red(:,:,:,round(end/2)),[],2))-minr));
            rgb(:,:,2) = squeeze(scg*(double(max(green(:,:,:,round(end/2)),[],2))-ming));
            rgb(:,:,3) = zeros(size(rgb(:,:,1)));
            if get(handles.takelog,'Value')
                rgb = log(double(rgb).*256+1)./log(256);
            end
            rgb2(:,:,1) = squeeze(scr*(double(max(red(:,:,:,round(end/2)),[],3))-minr));
            rgb2(:,:,2) = squeeze(scg*(double(max(green(:,:,:,round(end/2)),[],3))-ming));
            rgb2(:,:,3) = zeros(size(rgb2(:,:,1)));
            if get(handles.takelog,'Value')
                rgb2 = log(double(rgb2).*256+1)./log(256);
            end
        end
        % temp1 = squeeze(handles.inputdata(:,:,:,i));
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
        title('selectarea (depth - left, then right)');
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
    separateTiffs = get(handles.SeparateTiffs,'Value');
    skewbool = get(handles.Skewcorrection,'Value');
    delta = str2num(get(handles.SkewAng,'String'));
    m=0;
    conversionFactors = [handles.info.cal.ylat; handles.info.cal.zdep; handles.info.cal.xwid];
    if length(handles.ss)==3; handles.ss(4) = 1; end
    
    if handles.ss(4)==1
        if separateTiffs
            imgToSave1 = [FilePath,'uncorrected_RGB_' scanName, '_R.tiff'];
            imgToSave2 = [FilePath,'uncorrected_RGB_' scanName, '_G.tiff'];
            fileCounter = 1;
            while (2 == exist(imgToSave1, 'file'))
                imgToSave1 = [FilePath, 'uncorrected_RGB_' scanName '_R_' num2str(fileCounter) '.tiff'];
                imgToSave2 = [FilePath, 'uncorrected_RGB_' scanName '_G_' num2str(fileCounter) '.tiff'];
                fileCounter = fileCounter + 1;
            end
        else
            imgToSave = [FilePath,'uncorrected_RGB_' scanName, '.tiff'];
            fileCounter = 1;
            while (2 == exist(imgToSave, 'file'))
                imgToSave = [FilePath, 'RGB_' scanName '_' num2str(fileCounter) '.tiff'];
                fileCounter = fileCounter + 1;
            end
        end
        %normalize each volume to its max
        [red, green] = go_Callback(hObject, eventdata, handles,1);
        if crp0==0; crp0 = 1; crp1 = size(red,3); crp0y = 1; crp1y = size(red,1); crp2 = 1; crp3 = size(red,2); end
        red = squeeze(double(red(crp0y:crp1y,crp2:crp3,crp0:crp1)));
        green = squeeze(double(green(crp0y:crp1y,crp2:crp3,crp0:crp1)));
        
        if skewbool
            if separateTiffs
                imgToSave1 = [FilePath,'RGB_' scanName, '_R.tiff'];
                imgToSave2 = [FilePath,'RGB_' scanName, '_G.tiff'];
                fileCounter = 1;
                while (2 == exist(imgToSave1, 'file'))
                    imgToSave1 = [FilePath, 'RGB_' scanName '_R_' num2str(fileCounter) '.tiff'];
                    imgToSave2 = [FilePath, 'RGB_' scanName '_G_' num2str(fileCounter) '.tiff'];
                    fileCounter = fileCounter + 1;
                end
                
            else
                imgToSave = [FilePath,'RGB_' scanName, '.tiff'];
                fileCounter = 1;
                while (2 == exist(imgToSave, 'file'))
                    imgToSave = [FilePath, 'RGB_' scanName '_' num2str(fileCounter) '.tiff'];
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
            % Correct for skew
            %
            %             figure(99); subplot(2,1,1)
            %             imagesc([0:size(red, 1)-1]*conversionFactors(3),[0:size(red, 3)-1]*conversionFactors(2), log(squeeze(max(red(:,round(end/3):round(2*end/3), :), [], 2)))')
            %             colormap gray; axis image; title('Raw')
            
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
        
        
        
        if get(handles.ScaleTiffs,'Value')
            if get(handles.takelog,'Value')
                red = uint16(2^16*(log(scr*(double(red)-minr)+1)./log(2)));
                green = uint16(2^16*(log(scg*(double(green)-ming)+1)./log(2)));
                max(max(max(red)))
                min(min(min(red)))
            else
                red = uint16(2^16*squeeze(scr*(double(red)-minr)));
                green = uint16(2^16*squeeze(scg*(double(green)-ming)));
            end
        else
            if get(handles.takelog,'Value')
                red = uint16((log(double(red)+1)./log(2^14)));
                green = uint16((log(double(green)+1)./log(2^14)));
            else
                red = uint16(squeeze(red));
                green = uint16(squeeze(green));
            end
        end
        
        
        ss = size(red);
        for i = 1:ss(3)
            if separateTiffs
                imwrite(red(:,:,i), imgToSave1, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
                imwrite(green(:,:,i), imgToSave2, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
            else
                clear rgb
                rgb(:,:,1) = red(:,:,i);
                rgb(:,:,2) = green(:,:,i);
                rgb(:,:,3) = zeros(size(red(:,:,i)));
                temp11 = squeeze(rgb);%(2^16)*temp1/(max(max(max(temp1))));
                
                imwrite(temp11, imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
            end
        end
        
        if get(handles.MIP,'Value')
            tmpr = max(max(red(:,:,:),[],3)-min(red(:,:,:),[],3),[],2);
            [ar, br] = hist(reshape(double(tmpr),[1,prod(size(tmpr))]),200);
            cc1r = br(min(find(ar>0.01*mean(ar))));
            cc2r = br(max(find(ar>0.01*mean(ar))));
            tmpg = max(max(green(:,:,:),[],3)-min(green(:,:,:),[],3),[],2);
            [ag, bg] = hist(reshape(double(tmpg),[1,prod(size(tmpg))]),200);
            cc1g = bg(min(find(ag>0.01*mean(ag))));
            cc2g = bg(max(find(ag>0.01*mean(ag))));
            
            tempr = squeeze(max(red(:,2:end-2,:)-red(:,:,:),[],2))';
            tempr = cat(1,tempr,1000*ones(2,size(red,1)));
            tempr = cat(1,tempr, (squeeze(max(red(:,2:end-2,:)-red(:,:,:),[],3))'));
            tempr = uint8((256/(cc2r-cc1r))*(tempr-cc1r));
            tempg = squeeze(max(green(:,2:end-2,:)-green(:,:,:),[],2))';
            tempg = cat(1,tempg,1000*ones(2,size(green,1)));
            tempg = cat(1,tempg, (squeeze(max(green(:,2:end-2,:)-green(:,:,:),[],3))'));
            tempg = uint8((256/(cc2g-cc1g))*(tempg-cc1g));
            temprg = cat(3,tempr, tempg, zeros(size(tempr)));
            
            filename =  fullfile(dataDirectory, ['RGB' scanName '__MIP.jpg']);
            imwrite(tempy,filename,'jpg')
        end
    else
        %parfor_progress(handles.ss(4));
        %parfor (kk =1:handles.ss(4),8) %scans
        for kk =1:handles.ss(4)
            fileCounter = 1;
            if separateTiffs
                imgToSave1 = [FilePath,'uncorrected_RGB_' scanName '_t', num2str(kk) '_R.tiff'];
                imgToSave2 = [FilePath,'uncorrected_RGB_' scanName '_t', num2str(kk) '_G.tiff'];
                while (2 == exist(imgToSave1, 'file'))
                    imgToSave1 = [FilePath,'uncorrected_RGB_' scanName '_' num2str(fileCounter) '_t', num2str(kk) '_R.tiff'];
                    imgToSave2 = [FilePath,'uncorrected_RGB_' scanName '_' num2str(fileCounter) '_t', num2str(kk) '_G.tiff'];
                    fileCounter = fileCounter + 1;
                end
            else
                imgToSave = [FilePath,'uncorrected_RGB_' scanName, '_', num2str(kk) '.tiff'];
                while (2 == exist(imgToSave, 'file'))
                    imgToSave = [FilePath,'uncorrected_RGB_' scanName '_' num2str(fileCounter) '_t', num2str(kk) '.tiff'];
                    fileCounter = fileCounter + 1;
                end
            end
            
            %normalize each volume to its max
            [red, green] = go_Callback(hObject, eventdata, handles,kk);
            red = squeeze(double(red(crp0y:crp1y,crp2:crp3,crp0:crp1)));
            green = squeeze(double(green(crp0y:crp1y,crp2:crp3,crp0:crp1)));
            
            
            if skewbool
                if separateTiffs
                    imgToSave1 = [FilePath,'R_' scanName, '_', num2str(kk) '.tiff'];
                    imgToSave2 = [FilePath,'G_' scanName, '_', num2str(kk) '.tiff'];
                else
                    imgToSave = [FilePath,'RGB_' scanName, '_', num2str(kk) '.tiff'];
                end
                % Coordinate System Correction
                red = flip(flip(red(:,:,:),2),3);
                red = permute(red, [3 1 2]);
                RR = imref3d(size(red), conversionFactors(1), conversionFactors(3), conversionFactors(2));
                
                green = flip(flip(green(:,:,:),2),3);
                green = permute(green, [3 1 2]);
                RG = imref3d(size(green), conversionFactors(1), conversionFactors(3), conversionFactors(2));
                
                % Correct for skew
                %                     figure(99); subplot(1,2,1)
                %                     imagesc([0:size(red, 1)-1]*conversionFactors(3),[0:size(red, 3)-1]*conversionFactors(2), squeeze(max(red(:,round(end/3):round(2*end/3), :), [], 2))')
                %                     colormap gray; axis image
                
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
                         
         if get(handles.ScaleTiffs,'Value')
            if get(handles.takelog,'Value')
                red = uint16(2^16*(log(scr*(double(red)-minr+1))./log(2)));
                green = uint16(2^16*(log(scg*(double(green)-ming+1))./log(2)));
            else
                red = uint16(2^16*squeeze(scr*(double(red)-minr)));
                green = uint16(2^16*squeeze(scg*(double(green)-ming)));
            end
        else
            if get(handles.takelog,'Value')
                red = uint16((log(double(red)+1)./log(2^16)));
                green = uint16((log(double(green)+1)./log(2^16)));
            else
                red = uint16(squeeze(red));
                green = uint16(squeeze(green));
            end
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
                tmpr = max(max(red(:,:,:),[],3)-min(red(:,:,:),[],3),[],2);
                [ar, br] = hist(reshape(double(tmpr),[1,prod(size(tmpr))]),200);
                cc1r = br(min(find(ar>0.01*mean(ar))));
                cc2r = br(max(find(ar>0.01*mean(ar))));
                tmpg = max(max(green(:,:,:),[],3)-min(green(:,:,:),[],3),[],2);
                [ag, bg] = hist(reshape(double(tmpg),[1,prod(size(tmpg))]),200);
                cc1g = bg(min(find(ag>0.01*mean(ag))));
                cc2g = bg(max(find(ag>0.01*mean(ag))));
                
                tempr = squeeze(max(red(:,2:end-2,:)-red(:,:,:),[],2))';
                tempr = cat(1,tempr,1000*ones(2,size(red,1)));
                tempr = cat(1,tempr, (squeeze(max(red(:,2:end-2,:)-red(:,:,:),[],3))'));
                tempr = uint8((256/(cc2r-cc1r))*(tempr-cc1r));
                tempg = squeeze(max(green(:,2:end-2,:)-green(:,:,:),[],2))';
                tempg = cat(1,tempg,1000*ones(2,size(green,1)));
                tempg = cat(1,tempg, (squeeze(max(green(:,2:end-2,:)-green(:,:,:),[],3))'));
                tempg = uint8((256/(cc2g-cc1g))*(tempg-cc1g));
                temprg = cat(3,tempr, tempg, zeros(size(tempr)));
                M(kk) = im2frame(temprg);
            end
            
            %parfor_progress;
        end
        if get(handles.MIP, 'Value')
            writeVideo(vidObj,M);
            close(vidObj)
        end
        if skewbool
            %                     figure(99); subplot(1,2,2)
            %                     imagesc([0:size(red, 1)-1]*conversionFactors(3),[0:size(red, 3)-1]*conversionFactors(2), squeeze(max(red(:,round(end/3):round(2*end/3), :,1), [], 2))')
            %                     colormap gray; axis image
            
            disp('Skew correction applied')
        end
        
        % toc
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
handles.proc.scalered = scr;
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
set(handles.scalegreen2,'String',num2str(scg));
set(handles.autoscale,'Value',0);
handles.proc.scalegreen = scg;
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
    set(handles.scalegreen2,'String',handles.proc.autosc_scalegreen);
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
if (str2num(get(handles.RotateGreen,'String'))~=0)
    
    axis(handles.MainFigure);
    clear rgb
    rgb = handles.RGBMIPS;
    
    rotval = str2num(get(handles.RotateGreen,'String'))-handles.rotpastval;
    rgb(:,:,2) = imrotate(rgb(:,:,2),rotval,'bilinear','Crop');
    imagesc(squeeze(uint8((256)*rgb)));
    
    handles.rotpastval = str2num(get(handles.RotateGreen,'String'));
    
    guidata(hObject,handles)
    display('image adjusted')
    
end


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
[red, ~] = go_Callback(hObject, eventdata, handles,200);

% Coordinate System Correction
conversionFactors = [handles.info.cal.ylat; handles.info.cal.zdep; handles.info.cal.xwid];
red = flip(flip(red,2),3);
red = permute(red, [3 1 2]);
RR = imref3d(size(red), conversionFactors(1), conversionFactors(3), conversionFactors(2));

% Correct for skew
figure(99); subplot(2,1,1)
imagesc([0:size(red, 1)-1]*conversionFactors(3),[0:size(red, 3)-1]*conversionFactors(2), log(squeeze(double(max(red(:,round(end/3):round(2*end/3), :), [], 2))))')
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
imagesc([0:size(red, 1)-1]*conversionFactors(3),[0:size(red, 3)-1]*conversionFactors(2), log(double(squeeze(max(red(:,round(end/3):round(2*end/3), :), [], 2)))'))
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
transforms.scg = str2num(get(handles.scalegreen2,'String'));
transforms.minr = str2num(get(handles.minr,'String'));
transforms.ming = str2num(get(handles.ming,'String'));
transforms.xg = handles.proc.rgb.xg2;
transforms.xr =handles.proc.rgb.xr2;
transforms.yr = handles.proc.rgb.yr2;
transforms.yg = handles.proc.rgb.yg2;
transforms.zr = handles.proc.rgb.zr2;
transforms.zg = handles.proc.rgb.zg2;

FilePath = fullfile(handles.dataDirectory, 'tiff_stacks', 'RGBtransform.mat');

% checks for repeat file names to prevent overwrite of previous info files
fileCounter = 1;
while (2 == exist(FilePath, 'file'))
    FilePath =  fullfile(handles.dataDirectory, 'tiff_stacks', ['RGBtransform_' num2str(fileCounter) '.mat']);
    fileCounter = fileCounter + 1;
end

% save transform info file
save(FilePath, 'transforms')




% --- Executes on button press in LoadTransforms.
function LoadTransforms_Callback(hObject, eventdata, handles)
% hObject    handle to LoadTransforms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

presentDirectory = pwd;
cd(fullfile(handles.dataDirectory, 'tiff_stacks/'));
[fileName, pathName, ~] = uigetfile('*.mat');
cd (pathName)
transforms = load(fileName);
cd(presentDirectory)

set(handles.scalered,'String',num2str(transforms.transforms.scr));
set(handles.scalegreen2,'String',num2str(transforms.transforms.scg));
set(handles.minr,'String',num2str(transforms.transforms.minr));
set(handles.ming,'String',num2str(transforms.transforms.ming));
handles.proc.rgb.xg2 = transforms.transforms.xg;
handles.proc.rgb.xr2 = transforms.transforms.xr;
handles.proc.rgb.yr2 = transforms.transforms.yr;
handles.proc.rgb.yg2 = transforms.transforms.yg;
handles.proc.rgb.zr2 = transforms.transforms.zr;
handles.proc.rgb.zg2 = transforms.transforms.zg;
set(handles.RotateGreen,'String',num2str(transforms.transforms.rotateg));
set(handles.ScaleG,'String',num2str(transforms.transforms.scaleg));


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

% Hints: get(hObject,'String') returns contents of ScaleG as text
%        str2double(get(hObject,'String')) returns contents of ScaleG as a double

if (str2num(get(handles.ScaleG,'String'))~=0) || (str2num(get(handles.ScaleG,'String'))~=handles.scalepastval)
    
    scaleval = str2num(get(handles.ScaleG,'String'))/handles.scalepastval;
    
    if (str2num(get(handles.ScaleG,'String'))~=1)
        
        green = handles.RGBMIPS(:,:,2);
        green = imresize(green,scaleval);
        if str2num(get(handles.ScaleG,'String')) <1
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
    end
    
    handles.RGBMIPS(:,:,2) = green;
    
    axis(handles.MainFigure);
    imagesc(squeeze(uint8((256)*handles.RGBMIPS)));
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



function scalegreen2_Callback(hObject, eventdata, handles)
% hObject    handle to scalegreen2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scalegreen2 as text
%        str2double(get(hObject,'String')) returns contents of scalegreen2 as a double


% --- Executes during object creation, after setting all properties.
function scalegreen2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scalegreen2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in takelog.
function takelog_Callback(hObject, eventdata, handles)
% hObject    handle to takelog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of takelog
