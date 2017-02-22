function varargout = plotVol_v1(sig, n)
% PLOTVOL M-file for plotVol.fig
%      PLOTVOL, by itself, creates a new PLOTVOL or raises the existing
%      singleton*.
%
%      H = PLOTVOL returns the handle to a new PLOTVOL or the handle to
%      the existing singleton*.
%
%      PLOTVOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTVOL.M with the given input arguments.
%
%      PLOTVOL('Property','Value',...) creates a new PLOTVOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plotVol_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plotVol_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plotVol

% Last Modified by GUIDE v2.5 17-Feb-2011 15:57:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plotVol_OpeningFcn, ...
                   'gui_OutputFcn',  @plotVol_OutputFcn, ...
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


% --- Executes just before plotVol is made visible.
function plotVol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plotVol (see VARARGIN)

% Choose default command line output for plotVol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
if strcmp(get(hObject,'Visible'),'off')
    initialize_gui(handles);
    uiwait(hObject);
end
% UIWAIT makes plotVol wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plotVol_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% % % varargout{1} = handles.output;



function rot_x_Callback(hObject, eventdata, handles)
% hObject    handle to rot_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rot_x as text
%        str2double(get(hObject,'String')) returns contents of rot_x as a double
rx = str2double(get(hObject,'String'));
ry = str2double(get(handles.rot_y,'String'));
rz = str2double(get(handles.rot_z,'String'));
dValue = get(handles.slider1,'Value');
direccao = get(handles.listDir,'Value');
displayV(handles,rx,ry,rz,dValue,direccao);

% --- Executes during object creation, after setting all properties.
function rot_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rot_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rot_y_Callback(hObject, eventdata, handles)
% hObject    handle to rot_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rot_y as text
%        str2double(get(hObject,'String')) returns contents of rot_y as a double
ry = str2double(get(hObject,'String'));
rx = str2double(get(handles.rot_x,'String'));
rz = str2double(get(handles.rot_z,'String'));
dValue = get(handles.slider1,'Value');
direccao = get(handles.listDir,'Value');
displayV(handles,rx,ry,rz,dValue,direccao);

% --- Executes during object creation, after setting all properties.
function rot_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rot_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rot_z_Callback(hObject, eventdata, handles)
% hObject    handle to rot_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rot_z as text
%        str2double(get(hObject,'String')) returns contents of rot_z as a double
rz = str2double(get(hObject,'String'));
dValue = get(handles.slider1,'Value');
rx = str2double(get(handles.rot_x,'String'));
ry = str2double(get(handles.rot_y,'String'));
direccao = get(handles.listDir,'Value');
displayV(handles,rx,ry,rz,dValue,direccao);

% --- Executes during object creation, after setting all properties.
function rot_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rot_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
dValue = get(hObject,'Value'); 
rx = str2double(get(handles.rot_x,'String'));
ry = str2double(get(handles.rot_y,'String'));
rz = str2double(get(handles.rot_z,'String'));
direccao = get(handles.listDir,'Value');
displayV(handles,rx,ry,rz,dValue,direccao);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function initialize_gui(handles);


set(handles.slider1,'Max',1,'Min',-1,'Value',-1);
set(handles.listDir,'Value',3);

rx = str2double(get(handles.rot_x,'String'));
ry = str2double(get(handles.rot_y,'String'));
rz = str2double(get(handles.rot_z,'String'));
dValue = get(handles.slider1,'Value');
direccao = get(handles.listDir,'Value');
displayV(handles,rx,ry,rz,dValue,direccao);






% --- Executes on selection change in listDir.
function listDir_Callback(hObject, eventdata, handles)
% hObject    handle to listDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listDir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listDir
direccao = get(hObject,'Value');
rx = str2double(get(handles.rot_x,'String'));
ry = str2double(get(handles.rot_y,'String'));
rz = str2double(get(handles.rot_z,'String'));
dValue = get(handles.slider1,'Value');
displayV(handles,rx,ry,rz,dValue,direccao);


% --- Executes during object creation, after setting all properties.
function listDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function displayV(handles,rx,ry,rz,dValue,direccao)


%load temp.mat; % read sig and n

[x,y,z] = meshgrid(linspace(-1,1,n),linspace(-1,1,n),linspace(-1,1,n));
xdir = [1 0 0]; ydir = [0 1 0]; zdir = [0 0 1]; % rotate along these directions
center = [0 0 0];

axes(handles.vol),
% creates a default surface and rotates it

hsp = surf(linspace(-1,1,n),linspace(-1,1,n),zeros(n)+dValue);
switch direccao,
    case 1
        rotate(hsp,ydir,90,center)
    case 2
        rotate(hsp,xdir,90,center)
end

rotate(hsp,xdir,rx,center)
rotate(hsp,ydir,ry,center)
rotate(hsp,zdir,rz,center)

xd = get(hsp,'XData'); % plane coordinates
yd = get(hsp,'YData');
zd = get(hsp,'ZData');
delete(hsp)


% slice(x,y,z,sig,[-1,1],1,-1) % Draw some volume boundaries
% hh =  get(handles.vol,'Children');
% set(hh(1),'FaceAlpha',0.1);
% set(hh(2),'FaceAlpha',0.1);
% set(hh(3),'FaceAlpha',0.1);
% set(hh(4),'FaceAlpha',0.1);

cla
%  ----  Object surface representation  -----%
global v;
p=patch(isosurface(x,y,z,v,1));
set(p,'facecolor','red','edgecolor','none');
alpha(0.5); 
daspect([1,1,1])
hold on;
%------------------------%
h = slice(x,y,z,sig,xd,yd,zd);
set(h,'FaceAlpha',1);
title('Selected Slice');
pd = get(h,'CData'); % plane grey values

axes(handles.slice);
Title('Slice')
hold off
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]); axis square
view(-5,10),
axis off;

imshow(pd)
caxis([0 1])



