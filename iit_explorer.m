function varargout = iit_explorer(varargin)
% IIT_EXPLORER MATLAB code for iit_explorer.fig
%      IIT_EXPLORER, by itself, creates a new IIT_EXPLORER or raises the existing
%      singleton*.
%
%      H = IIT_EXPLORER returns the handle to a new IIT_EXPLORER or the handle to
%      the existing singleton*.
%
%      IIT_EXPLORER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IIT_EXPLORER.M with the given input arguments.
%
%      IIT_EXPLORER('Property','Value',...) creates a new IIT_EXPLORER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before iit_explorer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to iit_explorer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help iit_explorer

% Last Modified by GUIDE v2.5 12-Jul-2012 16:59:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @iit_explorer_OpeningFcn, ...
                   'gui_OutputFcn',  @iit_explorer_OutputFcn, ...
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


% --- Executes just before iit_explorer is made visible.
function iit_explorer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to iit_explorer (see VARARGIN)

% Choose default command line output for iit_explorer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes iit_explorer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = iit_explorer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in view_menu.
function view_menu_Callback(hObject, eventdata, handles)
% hObject    handle to view_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns view_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from view_menu


% --- Executes during object creation, after setting all properties.
function view_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nodes_list.
function nodes_list_Callback(hObject, eventdata, handles)
% hObject    handle to nodes_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nodes_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nodes_list


% --- Executes during object creation, after setting all properties.
function nodes_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nodes_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in complex_button.
function complex_button_Callback(hObject, eventdata, handles)
% hObject    handle to complex_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load fisheriris
plotmatrix(handles.overview_axes,meas)
linkdata on
brush on
