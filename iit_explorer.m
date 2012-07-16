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

% Last Modified by GUIDE v2.5 16-Jul-2012 11:33:13

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

% Set initial View
set(handles.Concepts,'Visible','Off')
set(handles.MIP,'Visible','Off')


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

view_choices = cellstr(get(hObject,'String'));
% turn all views off
for i = 1:length(view_choices)
    
    this_view = view_choices{i}(view_choices{i} ~= ' ');
    eval(['set(handles.' this_view ',''Visible'',''Off'')'])
    
end

% turn selected view on
selection = view_choices{get(hObject,'Value')};
selection = selection(selection ~= ' ');
eval(['set(handles.' selection ',''Visible'',''On'')'])


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
load('sample_partition_4n_sys.mat')
% conceptscatter(x,nWholeConcepts,handles.overview_axes,handles.mip_plot_panel);
[handles.mip_axes height extra_plots] = conceptscatter(x,nWholeConcepts,handles.overview_axes,handles.mip_plot_panel);
% setappdata(handles.mip_plot_panel,'PlotHeight',height,'ExtraPlots',extra_plots);
linkdata on



% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% set(handles.mip_plot_panel,'PlotHeight',height,'ExtraPlots',extra_plots);

% slide_val = get(hObject,'Value');
% % i could precombine these two values... we'll see...
% height = get(handles.mip_plot_panel,'PlotHeight');
% extra_plots = get(handles.mip_plot_panel,'ExtraPlots');
% offset = (height * extra_plots) * (1 - slide_val);
% for i= 1:length(handles.mip_axes)
%    
%     old_pos = get(ax{i},'Position');
% %     set(ax{i},'Position',old_pos + [0 offset 
%     
% end


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function brush_toggle_OffCallback(hObject, eventdata, handles)
% hObject    handle to brush_toggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
brush off


% --------------------------------------------------------------------
function brush_toggle_OnCallback(hObject, eventdata, handles)
% hObject    handle to brush_toggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
brush on

for i = 1:length(handles.mip_axes)
    
    brushed = findall(handles.mip_axes{i},'tag','Brushing');
    set(handles.brushed,'Parent',handles.mip_plot_panel,'Clipping','on')
    
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
