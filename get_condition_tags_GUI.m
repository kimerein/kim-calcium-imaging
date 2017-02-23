function varargout = get_condition_tags_GUI(varargin)
% GET_CONDITION_TAGS_GUI MATLAB code for get_condition_tags_GUI.fig
%      GET_CONDITION_TAGS_GUI, by itself, creates a new GET_CONDITION_TAGS_GUI or raises the existing
%      singleton*.
%
%      H = GET_CONDITION_TAGS_GUI returns the handle to a new GET_CONDITION_TAGS_GUI or the handle to
%      the existing singleton*.
%
%      GET_CONDITION_TAGS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GET_CONDITION_TAGS_GUI.M with the given input arguments.
%
%      GET_CONDITION_TAGS_GUI('Property','Value',...) creates a new GET_CONDITION_TAGS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before get_condition_tags_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to get_condition_tags_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help get_condition_tags_GUI

% Last Modified by GUIDE v2.5 22-Feb-2017 11:57:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @get_condition_tags_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @get_condition_tags_GUI_OutputFcn, ...
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


% --- Executes just before get_condition_tags_GUI is made visible.
function get_condition_tags_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to get_condition_tags_GUI (see VARARGIN)

% Choose default command line output for get_condition_tags_GUI
handles.output = hObject;

% Expected structures of these variables
% behForProfile = 
% 
%     [1x116 double]    [1x116 double]    [1x116 double]    [1x116 double]
%     [1x116 double]    [1x116 double]    [1x116 double]    [1x116 double]
%     [1x116 double]    [1x116 double]    [1x116 double]    [1x116 double]
%     [1x116 double]    [1x116 double]    [1x116 double]    [1x116 double]
%     [1x116 double]    [1x116 double]    [1x116 double]    [1x116 double]
%     [1x116 double]    [1x116 double]    [1x116 double]    [1x116 double]
%
% optoForProfile = 
% 
%     [  7x116 double]    [ 7x116 double]    [ 0x116 double]    [ 0x116 double]
%     [ 85x116 double]    [38x116 double]    [34x116 double]    [13x116 double]
%     [ 41x116 double]    [39x116 double]    [ 0x116 double]    [ 2x116 double]
%     [  1x116 double]    [ 1x116 double]    [ 0x116 double]    [ 0x116 double]
%     [ 92x116 double]    [45x116 double]    [34x116 double]    [13x116 double]
%     [134x116 double]    [85x116 double]    [34x116 double]    [15x116 double]
%
% out = 
% 
%     display_matrix: [24x46 double]
%           behavior: [1x1 struct]
%       optogenetics: [1x1 struct]
%             change: [1x1 struct]
%            sorting: [1x1 struct]
%               dist: [1x1 struct]
%             traces: [1x1 struct]
%           response: [1x1 struct]
%
% size(times)
% 
% ans =
% 
%      1   116
%
% withinCellAverages = 
% 
%     [46x116 double]    [46x116 double]    [46x116 double]    [46x116 double]
%     [46x116 double]    [46x116 double]    [46x116 double]    [46x116 double]
%     [46x116 double]    [46x116 double]    [46x116 double]    [46x116 double]
%     [46x116 double]    [46x116 double]    [46x116 double]    [46x116 double]
%     [46x116 double]    [46x116 double]    [46x116 double]    [46x116 double]
%     [46x116 double]    [46x116 double]    [46x116 double]    [46x116 double]
%
% withinCellResponses = 
% 
%     {1x46 cell}    {1x46 cell}    {1x46 cell}    {1x46 cell}
%     {1x46 cell}    {1x46 cell}    {1x46 cell}    {1x46 cell}
%     {1x46 cell}    {1x46 cell}    {1x46 cell}    {1x46 cell}
%     {1x46 cell}    {1x46 cell}    {1x46 cell}    {1x46 cell}
%     {1x46 cell}    {1x46 cell}    {1x46 cell}    {1x46 cell}
%     {1x46 cell}    {1x46 cell}    {1x46 cell}    {1x46 cell}
%
% withinCellStats = 
% 
%     [46x1 double]    [46x1 double]    [46x1 double]    [46x1 double]
%     [46x1 double]    [46x1 double]    [46x1 double]    [46x1 double]
%     [46x1 double]    [46x1 double]    [46x1 double]    [46x1 double]
%     [46x1 double]    [46x1 double]    [46x1 double]    [46x1 double]
%     [46x1 double]    [46x1 double]    [46x1 double]    [46x1 double]
%     [46x1 double]    [46x1 double]    [46x1 double]    [46x1 double]


in=varargin{1};
behForProfile=in.behForProfile;
optoForProfile=in.optoForProfile;
out=in.out;
times=in.times;
withinCellAverages=in.withinCellAverages;
withinCellResponses=in.withinCellResponses;
withinCellStats=in.withinCellStats;
optoMapping=in.optoMapping;
optoTimestamps=in.optoTimestamps;
optoStimTypes=in.optoStimTypes;
acq_obj=in.acq_obj;

handles.behForProfile=behForProfile;
handles.optoForProfile=optoForProfile;
handles.out=out;
handles.times=times;
handles.withinCellAverages=withinCellAverages;
handles.withinCellResponses=withinCellResponses;
handles.withinCellStats=withinCellStats;
handles.optoMapping=optoMapping;
handles.optoTimestamps=optoTimestamps;
handles.optoStimTypes=optoStimTypes;
handles.acq_obj=acq_obj;
handles.condition_labels=cell(size(withinCellAverages));

% Analysis code iterates on optoProfiles within behProfiles
% behProfiles are rows, optoProfiles are columns

set(hObject,'CloseRequestFcn',@closeGUI);

if size(optoMapping,1)~=size(withinCellAverages,2)
    error('optoMapping is not the same size as withinCellAverages');
end

k=1;
optoStimIntervals=[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5]; % in Volts, command to laser
for i=1:size(withinCellAverages,1)
    % Iterate through behavior profiles
    if strcmp(out.behavior.profile_type,'runningBeforeAndAfterOpto')
            behProfDescript=out.behavior.profiles{i};
            addToStr=[];
            switch i
                case 1
                    default_beh_string='norun_before_norun_after';
                case 2 
                    default_beh_string='run_before_run_after';
                case 3
                    default_beh_string='norun_before_run_after';
                case 4
                    default_beh_string='run_before_norun_after';
                case 5
                    default_beh_string='no_change_in_behavior';
                case 6
                    default_beh_string='all';
                case 7
                    default_beh_string='other';
            end
%             for l=1:length(behProfDescript)
%                 if isstr(behProfDescript{l})
%                     if strcmp(behProfDescript{l},'|')
%                         addToStr=[addToStr '_or'];
%                     elseif strcmp(behProfDescript{l},'&')
%                         addToStr=[addToStr '_and'];
%                     else
%                         addToStr=[addToStr '_unknownOperator'];
%                     end
%                 else
%                     curr=behProfDescript{l};
%                     if length(curr)~=2
%                         error('Unknown format for out.behavior -- look in analysisSettings.m for explanation of format');
%                     end
%                     if length(curr{1})==2
%                         addToStr=[addToStr 'both_before'];
%                     elseif curr{1}==0
%                         addToStr=[addToStr 'nonrun_before'];
%                     else
%                         addToStr=[addToStr 'run_before'];
%                     end
%                     if length(curr{2})==2
%                         addToStr=[addToStr 'both_after'];
%                     elseif curr{2}==0
%                         addToStr=[addToStr '_then_nonrun_after'];
%                     else
%                         addToStr=[addToStr '_then_run_after'];
%                     end
%                 end
%             end
%             default_beh_string=addToStr;
    else
        default_beh_string='';
    end
    for j=1:size(withinCellAverages,2)
        % Iterate through optogenetic profiles
        hax=subplot(size(withinCellAverages,1),size(withinCellAverages,2),k);
        axesname=['axes' num2str(k)];
        handles.(axesname)=hax;
        if strcmp(optoMapping{j,3},'control')
            plot(nanmean(optoTimestamps.timestamps_nocrop(ismember(optoStimTypes,j),:)- ...
                repmat(optoTimestamps.timestamps_nocrop(ismember(optoStimTypes,j),1),1,size(optoTimestamps.timestamps_nocrop,2)),1),...
                optoMapping{j,1},'Color','k');
        else
            plot(nanmean(optoTimestamps.timestamps_nocrop(ismember(optoStimTypes,j),:)- ...
                repmat(optoTimestamps.timestamps_nocrop(ismember(optoStimTypes,j),1),1,size(optoTimestamps.timestamps_nocrop,2)),1),...
                optoMapping{j,1},'Color','b');
        end
        hold on;
        % Plot behavior profile onto this same plot
        scaleBy=1000; % multiply behavior profile by this so similar scale to opto
        plot(times,behForProfile{i,j}*scaleBy,'Color',[0.8 0.7 0]);
        axis tight
        if j~=size(optoMapping,1)
            set(gca,'XTickLabelMode','manual','XTickLabel',[]);
        end
        set(gca,'YTick',[floor(max(optoMapping{j,1}))]);
        editname=['edit' num2str(j)];
        handles.(editname)=uicontrol('Style','edit','Parent',gcf,'Callback',@editCallback);
        [~,mi]=min(abs(max(optoMapping{j,1})/acq_obj.sabaMetadata.optoScaleFactor-optoStimIntervals));
        default_opto_string=['_Hz_' num2str(optoStimIntervals(mi)) 'V'];
        set(handles.(editname),'String',[default_beh_string default_opto_string '_' optoMapping{j,3}]); % set default name of this condition
        handles.condition_labels{i,j}=[default_beh_string default_opto_string '_' optoMapping{j,3}];
        set(gca,'units','pix');
        axesPosition=get(hax,'Position');
%         set(handles.(editname),'Position',[axesPosition(1)+axesPosition(3)+5 axesPosition(2) 200 20]);
        set(handles.(editname),'Position',[axesPosition(1)+axesPosition(3)-750 axesPosition(2)-30 750 20]);
        set(handles.(editname),'Tag',[num2str(i) '_' num2str(j)]);
        k=k+1;
    end
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes get_condition_tags_GUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = get_condition_tags_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
temp{1}=handles.condition_labels;
varargout{1} = temp;
delete(handles.figure1);

function  editCallback(hObject, eventdata)

[~,currFigure]=gcbo;
handles=guidata(currFigure);

currentString=get(hObject,'String');
whichEdit=get(hObject,'Tag');
underscoreInd=regexp(whichEdit,'_');
behNum=whichEdit(1:underscoreInd-1);
optoNum=whichEdit(underscoreInd+1:end);
handles.condition_labels{behNum,optoNum}=currentString;
disp(handles.condition_labels{behNum,optoNum});

% Update handles structure
guidata(hObject, handles);


function closeGUI(hObject, eventdata)

handles=guidata(hObject);

% % Update handles structure
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
