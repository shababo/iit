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

% Last Modified by GUIDE v2.5 19-Jul-2012 20:14:06

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

set(hObject, 'Renderer', 'painters')

% Choose default command line output for iit_explorer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes iit_explorer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Set initial View
set(handles.Concepts,'Visible','Off')
set(handles.MIP,'Visible','Off')

% load data
if nargin == 4 && isstruct(varargin{1})
    handles.data = varargin{1};
    guidata(hObject,handles)
else
    fprintf('Please load a struct when opening IIT_EXPLORER\n');
    handles.output = 'ERROR';
    guidata(hObject,handles)
    close(gcf)
    return
end

num_states = 2^handles.data.num_nodes;
all_states = ~isempty(handles.data.Big_phi_M{2});

% setup state listbox
if all_states
    states = cell(1,num_states + 1);
    states{1} = 'Average';
    for i = 1:num_states
        states{i + 1} = dec2bin(i-1,handles.data.num_nodes);
    end
else
    states = {mod_mat2str(handles.data.current_state')};
    set(handles.state_list,'Enable','off')
end
set(handles.state_list,'String',states);

% setup subset listbox
nodes = cell(1,handles.data.num_nodes);
for i = 1:handles.data.num_nodes
    nodes{i} = num2str(i);
end
set(handles.nodes_list,'String',nodes)
set(handles.nodes_list,'Value',handles.data.Complex{1})

set(handles.overview_axes_panel,'Visible','off');
set(handles.summary_panel,'Visible','off');



% set(handles.overview_scroll_panel,'Parent',handles.overview_axes_panel)





% --- Outputs from this function are returned to the command line.
function varargout = iit_explorer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


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
set(handles.mip_plot_panel,'Visible','off')
load('sample_partition.mat')
                                             
[handles.mip_axes height extra_plots] = conceptscatter(x,nWholeConcepts,handles.mip_main_axes,handles.mip_plot_panel);
% setappdata(handles.mip_plot_panel,'PlotHeight',height,'ExtraPlots',extra_plots);
linkdata on
set(handles.mip_plot_panel,'Visible','on')



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


% --- Executes on button press in refresh_button.
function refresh_button_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

view_choices = cellstr(get(handles.view_menu,'String'));
view = view_choices{get(handles.view_menu,'Value')};

subset = get(handles.nodes_list,'Value');
subset_index = convi(subset) - 1;

N = length(subset);

state_choice = get(handles.state_list,'Value');
if length(get(handles.state_list,'String')) == 1
    state_index = 1;
else
    state_index = trans10(flipud(trans2(state_choice - 2,handles.data.num_nodes)));
end

if strcmp(view,'Overview')
    
    set(handles.overview_axes_panel,'Visible','off')
    current_display_elements = allchild(handles.overview_axes_panel);

    for i = 1:length(current_display_elements)
            delete(current_display_elements(i))
    end



    set(handles.overview_axes_text,'Visible','off')
    set(handles.big_phi_text,'String',['Big Phi = ' num2str(handles.data.Big_phi_M{state_index}(subset_index))])
    set(handles.big_phi_MIP_text,'String',['Big Phi MIP = ' num2str(handles.data.Big_phi_MIP{state_index}(subset_index))])
    set(handles.MIP_text,'String',{'MIP:',[mod_mat2str(handles.data.complex_MIP_M{state_index}{subset_index}) '-'...
                                            mod_mat2str(pick_rest(subset,handles.data.complex_MIP_M{state_index}{subset_index}))]})
    set(handles.sum_small_phi_text,'String',['Sum Small Phi = ' num2str(sum(handles.data.small_phi_M{state_index}{subset_index}(:,1)))])
    set(handles.num_core_concepts_text,'String',['# Core Concepts = ' num2str(sum(handles.data.small_phi_M{state_index}{subset_index}(:,1) ~= 0))])

    [IRR_REP IRR_phi IRR_MIP M_IRR] = IRR_points(handles.data.concepts_M{state_index},...
                                                 handles.data.small_phi_M{state_index},...
                                                 handles.data.concept_MIP_M{state_index},subset, subset_index);


% 	disp(handles.overview_scroll_panel)
    plot_REP(handles.data.Big_phi_M{state_index}(subset_index), IRR_REP, IRR_phi, IRR_MIP,...
                                        handles.data.Complex{state_index}, handles.overview_axes_panel)

    % end    

    set(handles.summary_panel,'Visible','on')
    set(handles.overview_axes_panel,'Visible','on')
    set(handles.panel_slider,'Value',1.0)

elseif strcmp(view,'MIP')
   
    set(handles.mip_plot_panel,'Visible','off')
    set(handles.mip_loading_text,'Visible','on')
    current_display_elements = allchild(handles.mip_plot_panel);

    for i = 1:length(current_display_elements)
            delete(current_display_elements(i))
    end

    % get phi values for the whole
    w_phi_all = handles.data.small_phi_M{state_index}{subset_index}(:,1)';
    w_phi_concepts = w_phi_all(w_phi_all ~= 0);
%     IRR_whole = M_IRR_M{whole_i};
    
    
    % get concepts for the whole
    w_concept_dists_p = zeros(2^N,length(w_phi_concepts));
    w_concept_dists_f = zeros(2^N,length(w_phi_concepts));
    
    z = 1;
    for i = 1:length(w_phi_all)
        if (w_phi_all(i) ~= 0)
            
            if ~isempty(handles.data.concepts_M{state_index}{subset_index,1}{i}{1})
                w_concept_dists_p(:,z) = handles.data.concepts_M{state_index}{subset_index,1}{i}{1};
            end
            if ~isempty(handles.data.concepts_M{state_index}{subset_index,1}{i}{2})
                w_concept_dists_f(:,z) = handles.data.concepts_M{state_index}{subset_index,1}{i}{2};
            end
            z = z + 1;
        end
    end  
    
    
    
%     handles.data.concept_MIP_M{state_index}{subset_index} = concept_MIP_M_st;
    MIP_p1 = handles.data.complex_MIP_M{state_index}{subset_index};
    MIP_p1_index = convi(MIP_p1) - 1;
    MIP_p2 = pick_rest(subset,MIP_p1);
    MIP_p2_index = convi(MIP_p2) - 1;
    
    parts_phi_all = [handles.data.small_phi_M{state_index}{MIP_p1_index}(:,1)' ...
                          handles.data.small_phi_M{state_index}{MIP_p2_index}(:,1)'];
                                            
    nIRR = sum(parts_phi_all ~= 0);

    p_concept_dists_p = zeros(2^N,nIRR);
    p_concept_dists_f = zeros(2^N,nIRR);
    parts_phi_concepts = parts_phi_all(parts_phi_all ~= 0);

    z = 1;
    for k = 1:length(parts_phi_all)

        if (parts_phi_all(k) ~= 0)

            if(z <= sum(handles.data.small_phi_M{state_index}{MIP_p1_index}(:,1) ~= 0))
                p_concept_dists_p(:,z) = expand_prob(handles.data.concepts_M{state_index}{MIP_p1_index,1}{k}{1},subset,MIP_p1);
                p_concept_dists_f(:,z) = expand_prob(handles.data.concepts_M{state_index}{MIP_p1_index,1}{k}{2},subset,MIP_p1);
            else
                k_offset = k - size(handles.data.small_phi_M{state_index}{MIP_p1_index},1);
                p_concept_dists_p(:,z) = ...
                    expand_prob(handles.data.concepts_M{state_index}{MIP_p2_index,1}{k_offset}{1},subset,MIP_p2);
                p_concept_dists_f(:,z) = ...
                    expand_prob(handles.data.concepts_M{state_index}{MIP_p2_index,1}{k_offset}{2},subset,MIP_p2);
            end
            z = z + 1;

        end

    end
    
    all_concepts_p = [w_concept_dists_p'; p_concept_dists_p'];    
    [handles.mip_axes height extra_plots] = conceptscatter(all_concepts_p,size(w_concept_dists_p,2),handles.mip_plot_panel);

    linkdata on
    set(handles.mip_loading_text,'Visible','off')
    set(handles.mip_plot_panel,'Visible','on') 
    
    
end
    



% --- Executes on selection change in state_list.
function state_list_Callback(hObject, eventdata, handles)
% hObject    handle to state_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns state_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from state_list


% --- Executes during object creation, after setting all properties.
function state_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to state_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function panel_slider_Callback(hObject, eventdata, handles)
% hObject    handle to panel_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

shift = get(hObject,'Value')
position = get(handles.overview_axes_panel,'Position')
position(2) = (1 - position(4))*shift;
set(handles.overview_axes_panel,'Position',position)


% --- Executes during object creation, after setting all properties.
function panel_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to panel_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
