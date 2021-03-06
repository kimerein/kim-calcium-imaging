function varargout = groupOptoStimsGUI(varargin)
% GROUPOPTOSTIMSGUI MATLAB code for groupOptoStimsGUI.fig
%      GROUPOPTOSTIMSGUI, by itself, creates a new GROUPOPTOSTIMSGUI or raises the existing
%      singleton*.
%
%      H = GROUPOPTOSTIMSGUI returns the handle to a new GROUPOPTOSTIMSGUI or the handle to
%      the existing singleton*.
%
%      GROUPOPTOSTIMSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GROUPOPTOSTIMSGUI.M with the given input arguments.
%
%      GROUPOPTOSTIMSGUI('Property','Value',...) creates a new GROUPOPTOSTIMSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before groupOptoStimsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to groupOptoStimsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help groupOptoStimsGUI

% Last Modified by GUIDE v2.5 18-Jul-2016 09:23:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @groupOptoStimsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @groupOptoStimsGUI_OutputFcn, ...
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


% --- Executes just before groupOptoStimsGUI is made visible.
function groupOptoStimsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to groupOptoStimsGUI (see VARARGIN)

% Choose default command line output for groupOptoStimsGUI
handles.output = hObject;

optoMapping=varargin{1};
optoStimTypes=varargin{2};
times=varargin{3};
handles.optoMapping=optoMapping;
handles.optoStimTypes=optoStimTypes;

set(hObject,'CloseRequestFcn',@closeGUI);

beforeMapping=nan(size(optoMapping,1),1);
for i=1:size(optoMapping,1)
%     hax=subaxis(size(optoMapping,1),1,i,'Spacing',0.03,'Padding',0,'Margin',0);
    hax=subplot(size(optoMapping,1),1,i);
    axesname=['axes' num2str(i)];
    handles.(axesname)=hax;
    if strcmp(optoMapping{i,3},'control')
        plot(times,optoMapping{i,1},'Color','k');
    else
        plot(times,optoMapping{i,1},'Color','b');
    end
    axis tight
    if i~=size(optoMapping,1)
        set(gca,'XTickLabelMode','manual','XTickLabel',[]);
    end
    set(gca,'YTick',[floor(max(optoMapping{i,1}))]);
    editname=['edit' num2str(i)];
    handles.(editname)=uicontrol('Style','edit','Parent',gcf,'Callback',@editCallback);
    set(handles.(editname),'String',num2str(i));
    set(gca,'units','pix');
    axesPosition=get(hax,'Position');
    set(handles.(editname),'Position',[axesPosition(1)+axesPosition(3)+20 axesPosition(2) 20 20]);
    set(handles.(editname),'Tag',num2str(i));
    beforeMapping(i)=optoMapping{i,2};
end
handles.beforeMapping=beforeMapping;
handles.afterMapping=beforeMapping;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes groupOptoStimsGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = groupOptoStimsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.optoMapping;
varargout{2} = handles.optoStimTypes;
delete(handles.figure1);

function  editCallback(hObject, eventdata)

[~,currFigure]=gcbo;
handles=guidata(currFigure);

currentString=get(hObject,'String');
whichEdit=get(hObject,'Tag');
handles.afterMapping(str2num(whichEdit))=str2num(currentString);

% Update handles structure
guidata(hObject, handles);


function closeGUI(hObject, eventdata)

% When user closes GUI, fix optoStimTypes and remake optoMapping
handles=guidata(hObject);
beforeVals=handles.beforeMapping;
newVals=handles.afterMapping;
optoStimTypes=handles.optoStimTypes;
optoMapping=handles.optoMapping;

changeVals=zeros(length(beforeVals),length(optoStimTypes));
for i=1:length(beforeVals)
    currVal=beforeVals(i);
    changeVals(i,:)=optoStimTypes==currVal;
end

for i=1:size(changeVals,1)
    newVal=newVals(i);
    optoStimTypes(logical(changeVals(i,:)))=newVal;
end

% Remake optoMapping
diffVals=unique(newVals);
newOptoMapping=cell(length(diffVals),2);
for i=1:length(diffVals)
    combineVals=find(newVals==diffVals(i));
    sumOpto=zeros(1,length(optoMapping{1,1}));
    catName=[];
    for j=1:length(combineVals)
        sumOpto=sumOpto+optoMapping{combineVals(j),1};
        catName=strcat(catName,optoMapping{combineVals(j),3});
    end
    sumOpto=sumOpto./length(combineVals);
    newOptoMapping{i,1}=sumOpto;
    newOptoMapping{i,2}=diffVals(i);
    newOptoMapping{i,3}=catName;
end

handles.optoStimTypes=optoStimTypes;
handles.optoMapping=newOptoMapping;

% Update handles structure
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
