function varargout = iit(varargin)
% IIT MATLAB code for iit.fig
%      IIT, by itself, creates a new IIT or raises the existing
%      singleton*.
%
%      H = IIT returns the handle to a new IIT or the handle to
%      the existing singleton*.
%
%      IIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IIT.M with the given input arguments.
%
%      IIT('Property','Value',...) creates a new IIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before iit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to iit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help iit

% Last Modified by GUIDE v2.5 11-Jul-2012 15:54:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @iit_OpeningFcn, ...
                   'gui_OutputFcn',  @iit_OutputFcn, ...
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


% --- Executes just before iit is made visible.
function iit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to iit (see VARARGIN)

% Choose default command line output for iit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Set initial View
set(handles.Options,'Visible','Off')


% UIWAIT makes iit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = iit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function num_nodes_Callback(hObject, eventdata, handles)
% hObject    handle to num_nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_nodes as text
%        str2double(get(hObject,'String')) returns contents of num_nodes as a double

if isnan(str2double(get(hObject,'String'))) || ~isposintscalar(str2double(get(hObject,'String')))

    set(handles.warning,'String','Number of nodes must be a positive integer.');
    set(hObject,'String',num2str(size(get(handles.TPM,'Data'),1)));
    
else
    set(handles.warning,'String','');
    tpm_choices = cellstr(get(handles.tpm_type_menu,'String'));
    tpm_choice = tpm_choices{get(handles.tpm_type_menu,'Value')};
    updateTPMview(handles, tpm_choice)
    updateCurrentStateView(handles)
end


% --- Executes during object creation, after setting all properties.
function num_nodes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in net_definition_method.
function net_definition_method_Callback(hObject, eventdata, handles)
% hObject    handle to net_definition_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns net_definition_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from net_definition_method


% --- Executes during object creation, after setting all properties.
function net_definition_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to net_definition_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% redraw the TPM based on current settings
function updateTPMview(handles, tpm_choice)

nNodes = str2double(get(handles.num_nodes,'String'));

tpm_old = get(handles.TPM,'Data');
tpm_size_old = size(tpm_old,2);

if strcmp(tpm_choice,'State X State')
    
    % update TPM

    tpm_size_new = 2^nNodes;
    % increase size
    if tpm_size_new > tpm_size_old

        tpm_new = eye(tpm_size_new);
        tpm_new(1:size(tpm_old,1),1:size(tpm_old,2)) = tpm_old;

    % decrease size
    else

        % resize
        tpm_new = tpm_old(1:tpm_size_new,1:tpm_size_new);

    end

    set(handles.TPM,'Data',tpm_new);

    % rename cols and rows
    names = cell(1,tpm_size_new);
    for i = 1:tpm_size_new
        names{i} = dec2bin(i-1,nNodes);
    end
    set(handles.TPM,'ColumnName',names,'RowName',names);
    set(handles.TPM,'ColumnEditable',true(1,tpm_size_new));

elseif strcmp(tpm_choice,'State X Node')

    
    if nNodes > tpm_size_old
        
        tpm_new = zeros(2^nNodes,nNodes);
        tpm_new(1:2^tpm_size_old,1:tpm_size_old) = tpm_old;
    else
        tpm_new = tpm_old(1:2^nNodes,1:nNodes);
    end
    
    set(handles.TPM,'Data',tpm_new);
    
    row_names = cell(1,2^nNodes);
    for i = 1:2^nNodes
        row_names{i} = dec2bin(i-1,nNodes);
    end
    set(handles.TPM,'RowName',row_names);
    col_names = cell(1,nNodes);
    for i = 1:nNodes
        col_names{i} = num2str(i);
    end
    set(handles.TPM,'ColumnName',col_names);
    set(handles.TPM,'ColumnEditable',true(1,nNodes));
end



    
    

% redraw the current state based on new input
function updateCurrentStateView(handles)

nNodes = str2double(get(handles.num_nodes,'String'));

% update current state table

cur_state_old = get(handles.cur_state,'Data');
cur_state_size_old = length(cur_state_old);

% increase size
if nNodes > cur_state_size_old
    
    cur_state_new = zeros(1,nNodes);
    cur_state_new(1:cur_state_size_old) = cur_state_old;

% decrease size
else
    
    cur_state_new = cur_state_old(1:nNodes);
    
end

set(handles.cur_state,'Data',cur_state_new);
set(handles.cur_state,'ColumnEditable',true(1,nNodes));



% --- Executes on selection change in view_select.
function view_select_Callback(hObject, eventdata, handles)
% hObject    handle to view_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns view_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from view_select

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
function view_select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in parallel_option_menu.
function parallel_option_menu_Callback(hObject, eventdata, handles)
% hObject    handle to parallel_option_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns parallel_option_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parallel_option_menu


% --- Executes during object creation, after setting all properties.
function parallel_option_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parallel_option_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in state_option_menu.
function state_option_menu_Callback(hObject, eventdata, handles)
% hObject    handle to state_option_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns state_option_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from state_option_menu


% --- Executes during object creation, after setting all properties.
function state_option_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to state_option_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in complex_option_menu.
function complex_option_menu_Callback(hObject, eventdata, handles)
% hObject    handle to complex_option_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns complex_option_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from complex_option_menu


% --- Executes during object creation, after setting all properties.
function complex_option_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to complex_option_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% THIS WILL GET CHANGED WHEN OPTIONS ARE NAILED DOWN

op_big_phi = get(handles.big_phi_alg_menu,'Value');
if op_big_phi < 3
    op_big_phi = op_big_phi - 1;
end

op_normalize = get(handles.normalization_menu,'Value');
op_normalize_small_phi = 0;
op_normalize_big_phi = 0;
if op_normalize == 1 || op_normalize == 2
    op_normalize_small_phi = 1;
end
if op_normalize == 1 || op_normalize == 3
    op_normalize_big_phi = 1;
end

op_complex = get(handles.complex_option_menu,'Value');
if op_complex == 2
    op_complex = 0;
end

op_small_phi = get(handles.small_phi_func_menu,'Value') - 1;
op_big_phi_dist = get(handles.big_phi_func_menu,'Value') - 1;

op_ave = get(handles.state_option_menu,'Value') - 1;
op_parallel = get(handles.parallel_option_menu,'Value') - 1;

options = [3 1 2 1 1 0 0 1 1 1 op_big_phi 0 ...
           op_normalize_big_phi op_normalize_small_phi op_complex op_small_phi op_big_phi_dist op_ave op_parallel];
       

tpm_choices = cellstr(get(handles.tpm_type_menu,'String'));
tpm_choice = tpm_choices{get(handles.tpm_type_menu,'Value')};       

tpm = get(handles.TPM,'Data');

if strcmp(tpm_choice,'State X State')
    num_states = size(tpm,1);
    num_nodes = str2double(get(handles.num_nodes,'String'));
    new_tpm = zeros(num_states,num_nodes);
    
    for i = 1:num_states
        for j = 1:num_nodes
            for k = 1:num_states
                
                state = dec2bin(k-1,num_nodes);
                
                if strcmp(state(j),'1')
                    new_tpm(i,j) = new_tpm(i,j) + tpm(i,k);
                end
            end
        end
    end
    
    tpm = new_tpm;
end

current_state = get(handles.cur_state,'Data')';
noise = str2double(get(handles.noise,'String'));

iit_run(tpm,current_state,noise,options);



% --- Executes on selection change in big_phi_alg_menu.
function big_phi_alg_menu_Callback(hObject, eventdata, handles)
% hObject    handle to big_phi_alg_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns big_phi_alg_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from big_phi_alg_menu


% --- Executes during object creation, after setting all properties.
function big_phi_alg_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to big_phi_alg_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in big_phi_func_menu.
function big_phi_func_menu_Callback(hObject, eventdata, handles)
% hObject    handle to big_phi_func_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns big_phi_func_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from big_phi_func_menu


% --- Executes during object creation, after setting all properties.
function big_phi_func_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to big_phi_func_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in small_phi_func_menu.
function small_phi_func_menu_Callback(hObject, eventdata, handles)
% hObject    handle to small_phi_func_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns small_phi_func_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from small_phi_func_menu


% --- Executes during object creation, after setting all properties.
function small_phi_func_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to small_phi_func_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noise_Callback(hObject, eventdata, handles)
% hObject    handle to noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noise as text
%        str2double(get(hObject,'String')) returns contents of noise as a double

if isnan(str2double(get(hObject,'String')))
    set(handles.warning,'String','Noise must be a real number in [0,.5].');
    set(hObject,'String','0');
else
    noise = str2double(get(hObject,'String'));
    if noise < 0
        set(handles.warning,'String','Noise must be in [0,.5].');
        set(hObject,'String','0');
    elseif noise > .5
        set(handles.warning,'String','Noise must be in [0,.5].');
        set(hObject,'String','.5');
    else
        set(handles.warning,'String','');
    end
end

% --- Executes during object creation, after setting all properties.
function noise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in normalization_menu.
function normalization_menu_Callback(hObject, eventdata, handles)
% hObject    handle to normalization_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns normalization_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from normalization_menu


% --- Executes during object creation, after setting all properties.
function normalization_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to normalization_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in TPM.
function TPM_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to TPM (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

tpm = get(hObject,'Data');
tpm_choices = cellstr(get(handles.tpm_type_menu,'String'));
tpm_choice = tpm_choices{get(handles.tpm_type_menu,'Value')};

if any(isnan(eventdata.NewData) || any(eventdata.NewData < 0) || any(eventdata.NewData > 1))
    
    set(handles.warning,'String','Entries in the TPM must be real numbers in [0,1].');
    tpm(eventdata.Indices(1),eventdata.Indices(2)) = eventdata.PreviousData;
    set(hObject,'Data',tpm);
    
elseif strcmp(tpm_choice,'State X State') && any(sum(tpm,2) ~= 1)
    
    set(handles.warning,'String','Rows in the TPM must sum to 1');
    
else
    
    set(handles.warning,'String','');
    
end

    


% --- Executes when entered data in editable cell(s) in cur_state.
function cur_state_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to cur_state (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

state_vec = get(hObject,'Data');
if (eventdata.NewData ~= 0 && eventdata.NewData ~= 1)
    
    set(handles.warning,'String','Node states can only be 0 or 1')
    state_vec(eventdata.Indices(2)) = eventdata.PreviousData;
    set(hObject,'Data',state_vec)
else
    
    set(handles.warning,'String','');
    
end


% --- Executes on button press in upload_tpm.
function upload_tpm_Callback(hObject, eventdata, handles)
% hObject    handle to upload_tpm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = uigetfile('*.mat');

if (filename ~= 0)
    
    load(filename)

    if exist('tpm')
        if size(tpm,2) == size(tpm,1) && isposintscalar(log2(size(tpm,1)))

            num_nodes = log2(size(tpm,1));

            set(handles.num_nodes,'String',num2str(num_nodes))
            set(handles.tpm_type_menu,'Value',1); % state x state
set(handles.TPM,'Data',tpm);
            updateTPMview(handles,'State X State');
            

        elseif size(tpm,1) == 2^size(tpm,2)

            num_nodes = size(tpm,2);

            set(handles.num_nodes,'String',num2str(num_nodes))
            set(handles.tpm_type_menu,'Value',2); % state x node
            set(handles.TPM,'Data',tpm);

            updateTPMview(handles,'State X Node');
        end

        updateCurrentStateView(handles)
        
        set(handles.warning,'String','');
    else
        set(handles.warning,'String','No variable named ''tpm'' in that data file.')
    end
    
end
    

% --- Executes on selection change in tpm_type_menu.
function tpm_type_menu_Callback(hObject, eventdata, handles)
% hObject    handle to tpm_type_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tpm_type_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tpm_type_menu

tpm_choices = cellstr(get(hObject,'String'));
tpm_choice = tpm_choices{get(hObject,'Value')};
tpm = get(handles.TPM,'Data');
num_nodes = str2double(get(handles.num_nodes,'String'));
num_states = 2^num_nodes;

% if we want state x state and we were in state x node
if strcmp(tpm_choice,'State X State') && size(tpm,2) == num_nodes
    
    new_tpm = ones(num_states);
    
    % create state x state tpm
    for i = 1:num_states
        for j = 1:num_states
            
            state = num2str(dec2bin(j-1,num_nodes));
            
            for k = 1:num_nodes
                if strcmp(state(k),'1')
                    new_tpm(i,j) = new_tpm(i,j) * tpm(i,k);
                else
                    new_tpm(i,j) = new_tpm(i,j) * (1 - tpm(i,k));
                end
            end
        end
    end
% if we want state x node and we were in state x state    
elseif strcmp(tpm_choice,'State X Node') && size(tpm,2) == num_states
    
    new_tpm = zeros(num_states,num_nodes);
    
    for i = 1:num_states
        for j = 1:num_nodes
            for k = 1:num_states
                
                state = dec2bin(k-1,num_nodes);
                
                if strcmp(state(j),'1')
                    new_tpm(i,j) = new_tpm(i,j) + tpm(i,k);
                end
            end
        end
    end
end

updateTPMview(handles, tpm_choice)
set(handles.TPM,'Data',new_tpm)


% --- Executes during object creation, after setting all properties.
function tpm_type_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tpm_type_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function NetworkDefinition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NetworkDefinition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'Max',2)
