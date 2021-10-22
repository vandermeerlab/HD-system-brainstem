function varargout = seqxplore(varargin)
% SEQXPLORE M-file for seqxplore.fig
%      SEQXPLORE, by itself, creates a new SEQXPLORE or raises the existing
%      singleton*.
%
%      H = SEQXPLORE returns the handle to a new SEQXPLORE or the handle to
%      the existing singleton*.
%
%      SEQXPLORE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEQXPLORE.M with the given input arguments.
%
%      SEQXPLORE('Property','Value',...) creates a new SEQXPLORE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before seqxplore_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to seqxplore_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help seqxplore

% Last Modified by GUIDE v2.5 21-Jun-2011 17:38:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @seqxplore_OpeningFcn, ...
                   'gui_OutputFcn',  @seqxplore_OutputFcn, ...
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


% --- Executes just before seqxplore is made visible.
function seqxplore_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to seqxplore (see VARARGIN)

% Choose default command line output for seqxplore
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes seqxplore wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = seqxplore_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadButton.
function LoadButton_Callback(hObject, eventdata, handles)
% pop up dialog for file loading, draw permanents (pos samples, markers,
% rat)

% load
[input_file,pathname] = uigetfile( ...
       {'*.mat', 'MATLAB Data Files (*.mat)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Select files', ... 
        'MultiSelect', 'off');
 
if pathname == 0
    return
end

load(fullfile(pathname,input_file));
handles.pf = pf; clear pf;
if exist('cp','var') % recon data
    handles.cp = cp; clear cp;
end

whitebg(gcf,'k');

% plot bg
plotPosData(Data(handles.pf.sd_orig.x),Data(handles.pf.sd_orig.y),handles.TopLeftAxes);
plotPosData(Data(handles.pf.sd_orig.x),Data(handles.pf.sd_orig.y),handles.TopRightAxes);

% find start of first lap, put this time in time field
handles.xR = Range(handles.pf.sd_orig.x);
handles.lap_ids = nan(size(handles.xR));
for ixR = 1:length(handles.xR), handles.lap_ids(ixR) = sum(handles.pf.lapData.startTime <= handles.xR(ixR)); end;
firstLap = sort(find(handles.lap_ids == 1),'ascend'); firstLap = handles.xR(firstLap(1));

set(handles.CurrentEdit,'String',num2str(firstLap));

handles.xD = Data(handles.pf.sd_orig.x);
handles.yD = Data(handles.pf.sd_orig.y);

% set noCP
set(handles.CPmax,'String',num2str(length(handles.pf.cp_entries)));

% set noSW
set(handles.SWmax,'String',num2str(length(handles.pf.swt)));

% output switch cp
if isfield(handles.pf.sd_orig.ExpKeys,'SwitchT')
    [val,ind] = min(abs(handles.pf.cp_entries-handles.pf.sd_orig.ExpKeys.SwitchT));
    set(handles.CPmax,'String',cat(2,get(handles.CPmax,'String'),'(',num2str(ind),')'));
end

% initialize rat handle for MazeView
t = str2double(get(handles.CurrentEdit,'String')); % this is now a time
% need time to xR/xD conversion function
t = TimeToInd(t,handles.xR);

switch handles.pf.lapData.bottomChoice{handles.lap_ids(t)} % get lap type for current sample
    case 'L'
        axes(handles.TopLeftAxes); hold on;
        handles.rat = plot(handles.xD(t),handles.yD(t),'wo','MarkerSize',15);
    case 'R'
        axes(handles.TopRightAxes); hold on;
        handles.rat = plot(handles.xD(t),handles.yD(t),'wo','MarkerSize',15);
end

% draw some cells
% first, go through fields to get cell id and rank
field_count_L = 0;
field_count_R = 0;

for iC = 1:length(handles.pf.l.fields)
    temp = handles.pf.l.fields{iC};
    while ~isempty(temp)
       
        field_count_L = field_count_L + 1;
        handles.field_loc_L(field_count_L) = temp(1);
        handles.field_cell_L(field_count_L) = iC;
        
        if length(temp) > 1
        temp = temp(2:end);
        else
            temp = [];
        end
        
    end
    
    temp = handles.pf.r.fields{iC};
    while ~isempty(temp)
       
        field_count_R = field_count_R + 1;
        handles.field_loc_R(field_count_R) = temp(1);
        handles.field_cell_R(field_count_R) = iC;
        
        if length(temp) > 1
        temp = temp(2:end);
        else
            temp = [];
        end
        
    end
end

[out,tempL] = sort(handles.field_loc_L,'ascend');
[out,tempR] = sort(handles.field_loc_R,'ascend');

for iF = 1:length(tempL)
    handles.field_rank_L(tempL(iF)) = iF;
end

for iF = 1:length(tempR)
    handles.field_rank_R(tempR(iF)) = iF;
end

% then, make some colors -- fields from the same cell should have same
% color
nCols = max(cat(2,handles.field_cell_L,handles.field_cell_R));
temp_cols = cat(2,rand(nCols,1),rand(nCols,1),rand(nCols,1));
for iF = 1:length(handles.field_rank_L)
    handles.field_col_L(iF,:) = temp_cols(handles.field_cell_L(iF),:);
end
for iF = 1:length(handles.field_rank_R)
    handles.field_col_R(iF,:) = temp_cols(handles.field_cell_R(iF),:);
end

% make tvec that gives which cp/ct pass we're in
handles.cpp = nan(size(handles.xR));
for iCP = 1:length(handles.pf.cp_entries)
    goodInd = find(handles.xR > handles.pf.cp_entries(iCP)-1 & handles.xR < handles.pf.cp_exits(iCP)+2);
    handles.cpp(goodInd) = iCP;
end

% plot CSCs
handles = plotLFPs(hObject,handles);

% now actually plot something
t = str2double(get(handles.CurrentEdit,'String'));
handles = plotMazeCells(hObject,handles,t);
handles = plotFieldCells(hObject,handles);
handles = centerFieldCells(hObject,handles,t);




% aaand recon
handles = plotRecon(hObject,handles,t);

% Update handles structure
guidata(hObject, handles);

function StepSizeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to StepSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StepSizeEdit as text
%        str2double(get(hObject,'String')) returns contents of StepSizeEdit as a double


% --- Executes during object creation, after setting all properties.
function StepSizeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StepSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CurrentEdit_Callback(hObject, eventdata, handles)
% hObject    handle to CurrentEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrentEdit as text
%        str2double(get(hObject,'String')) returns contents of CurrentEdit as a double


% --- Executes during object creation, after setting all properties.
function CurrentEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrentEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GoButton.
function GoButton_Callback(hObject, eventdata, handles)
% hObject    handle to GoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% find out what sample we are on
t = str2double(get(handles.CurrentEdit,'String'));

% update rat location
handles = updateRat(hObject,handles,t);
handles = plotMazeCells(hObject,handles,t);
handles = centerFieldCells(hObject,handles,t);
handles = plotRecon(hObject,handles,t);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in FastRevButton.
function FastRevButton_Callback(hObject, eventdata, handles)
% hObject    handle to FastRevButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SlowRevButton.
function SlowRevButton_Callback(hObject, eventdata, handles)
% hObject    handle to SlowRevButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SlowFwdButton.
function SlowFwdButton_Callback(hObject, eventdata, handles)
% hObject    handle to SlowFwdButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% find out what sample we are on
t = str2double(get(handles.CurrentEdit,'String'));
dt = str2double(get(handles.StepSizeEdit,'String'));
t_new = t+dt;
set(handles.CurrentEdit,'String',num2str(t_new));

% update rat location
handles = updateRat(hObject,handles,t_new);
handles = plotMazeCells(hObject,handles,t_new);
handles = centerFieldCells(hObject,handles,t_new);
handles = plotRecon(hObject,handles,t);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in FastFwdButton.
function FastFwdButton_Callback(hObject, eventdata, handles)
% hObject    handle to FastFwdButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% find out what sample we are on
t = str2double(get(handles.CurrentEdit,'String'));
dt = str2double(get(handles.StepSizeEdit,'String'));
t_new = t+10*dt;
set(handles.CurrentEdit,'String',num2str(t_new));

% update rat location
handles = updateRat(hObject,handles,t_new);
handles = plotMazeCells(hObject,handles,t_new);
handles = centerFieldCells(hObject,handles,t_new);
handles = plotRecon(hObject,handles,t);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in PlayButton.
function PlayButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlayButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in StopButton.
function StopButton_Callback(hObject, eventdata, handles)
% hObject    handle to StopButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Plots background position sample data for MazeView
function plotPosData(x,y,axesName)

cla(axesName); %clear the axes
axes(axesName); %set the axes to plot
plot(x,y,'.','MarkerSize',1,'Color',[0.5 0.5 0.5]);
axis off;

function handles = updateRat(hObject,handles,t)

t = TimeToInd(t,handles.xR);
switch handles.pf.lapData.bottomChoice{handles.lap_ids(t)} % get lap type for current sample
    case 'L'
        set(handles.rat,'Parent',handles.TopLeftAxes,'XData',handles.xD(t),'YData',handles.yD(t));
    case 'R'
        set(handles.rat,'Parent',handles.TopRightAxes,'XData',handles.xD(t),'YData',handles.yD(t));
end

% Update handles structure
guidata(hObject, handles);

function handles = plotMazeCells(hObject,handles,t)

if isfield(handles,'cellLeftH')
    if ~isempty(handles.cellLeftH)
        delete(handles.cellLeftH);
    end

    if ~isempty(handles.cellRightH)
        delete(handles.cellRightH);
    end
end

% go through fields and get #spikes for each. if >0, plot (at field peak location) and create
% handle
nHandle = 0;
handles.cellLeftH = [];
dSample = 0.1; % sec each way
t = TimeToInd(t,handles.xR);

axes(handles.TopLeftAxes); hold on;
for iF = 1:length(handles.field_loc_L)
   
    spks = Data(Restrict(handles.pf.sd_orig.S{handles.field_cell_L(iF)},handles.xR(t)-dSample,handles.xR(t)+dSample));
    if ~isempty(spks)
        nHandle = nHandle + 1;
        zp = repmat(handles.field_loc_L(iF),[length(spks) 1]);
        if isempty(zp)
            continue;
        end
        coords = handles.pf.l.binToCoordL(:,zp);
        coords = coords+randn(size(coords));
        handles.cellLeftH(nHandle) = plot(coords(1,:),coords(2,:),'.','Color',handles.field_col_L(iF,:));
    end
    
end

nHandle = 0;
handles.cellRightH = [];
axes(handles.TopRightAxes); hold on;
for iF = 1:length(handles.field_loc_R)
   
    spks = Data(Restrict(handles.pf.sd_orig.S{handles.field_cell_R(iF)},handles.xR(t)-dSample,handles.xR(t)+dSample));
    if ~isempty(spks)
        nHandle = nHandle + 1;
        zp = repmat(handles.field_loc_R(iF),[length(spks) 1]);
        coords = handles.pf.r.binToCoordR(:,zp);
        coords = coords+randn(size(coords));
        handles.cellRightH(nHandle) = plot(coords(1,:),coords(2,:),'.','Color',handles.field_col_R(iF,:));
    end
    
end

% Update handles structure
guidata(hObject, handles);

function handles = plotLFPs(hObject,handles)

% plot HC fissure lfp in striatum-left
nCells = length(handles.pf.sd_str.S);

if iscell(handles.pf.sd_orig.ExpKeys.goodTheta)
    fc = handles.pf.sd_orig.ExpKeys.goodTheta{1};
else
    fc = handles.pf.sd_orig.ExpKeys.goodTheta;
end
csc = LoadCSC(fc,'ConvertToS',1);
csc = Restrict(csc,handles.pf.sd_orig.ExpKeys.TimeOnTrack,handles.pf.sd_orig.ExpKeys.TimeOffTrack);
cscR = Range(csc); cscD = Data(csc); clear csc;
cscD = decimate(cscD,10); cscR = cscR(1:10:end);
Fs = 1./median(diff(cscR));
cscD = locdetrend(cscD,Fs,[1 0.5]);

fb = [5 15]; % filter in (roughly) theta range
[b,a] = cheby1(4,0.5,[fb(1)/(Fs/2) fb(2)/(Fs/2)]);
cscD = filtfilt(b,a,cscD);

% rescale
cscD = cscD./(max(abs(cat(1,max(cscD),min(cscD)))));
cscD = cscD.*(nCells/2)+nCells/2;

axes(handles.LeftStr); hold on;
plot(cscR,cscD,'Color',[0 0.3 0.3]);

% plot VS lfp in striatum-right
nCells = length(handles.pf.sd_str.S);

if iscell(handles.pf.sd_orig.ExpKeys.goodGamma)
    fc = handles.pf.sd_orig.ExpKeys.goodGamma{1};
else
    fc = handles.pf.sd_orig.ExpKeys.goodGamma;
end
csc = LoadCSC(fc,'ConvertToS',1);
csc = Restrict(csc,handles.pf.sd_orig.ExpKeys.TimeOnTrack,handles.pf.sd_orig.ExpKeys.TimeOffTrack);
cscR = Range(csc); cscD = Data(csc); clear csc;
cscD = decimate(cscD,9); cscR = cscR(1:9:end);
Fs = 1./median(diff(cscR));
cscD = locdetrend(cscD,Fs,[1 0.5]);

fb = [5 100]; % filter in (roughly) theta range
[b,a] = cheby1(4,0.5,[fb(1)/(Fs/2) fb(2)/(Fs/2)]);
cscD = filtfilt(b,a,cscD);

% rescale
cscD = cscD./(max(abs(cat(1,max(cscD),min(cscD)))));
cscD = cscD.*(nCells/2)+nCells/2;

axes(handles.RightStr); hold on;
plot(cscR,cscD,'Color',[0.3 0 0.3]);

guidata(hObject, handles);

function handles = plotFieldCells(hObject,handles)

% first, plot striatum (left panel)
axes(handles.LeftStr); hold on; axis off;
nCells = length(handles.pf.sd_str.S);
for iC = 1:nCells
   spks = Data(handles.pf.sd_str.S{iC});
   plot(spks,iC.*ones(size(spks)),'.','Color',handles.field_col_L(end-iC,:));     
end

% set vertical lines at current location and window edges
handles.curLineSl = line([handles.xR(1) handles.xR(1)],[1 nCells]); set(handles.curLineSl,'Color',[0.5 0.5 0.5]);
handles.fwdLineSl = line([handles.xR(1)+0.1 handles.xR(1)+0.1],[1 nCells]); set(handles.fwdLineSl,'Color',[0.3 0.3 0.3]);
handles.revLineSl = line([handles.xR(1)-0.1 handles.xR(1)-0.1],[1 nCells]); set(handles.revLineSl,'Color',[0.3 0.3 0.3]);

% first, plot striatum (right panel)
axes(handles.RightStr); hold on; axis off;
nCells = length(handles.pf.sd_str.S);
for iC = 1:nCells
   spks = Data(handles.pf.sd_str.S{iC});
   plot(spks,iC.*ones(size(spks)),'.','Color',handles.field_col_L(end-iC,:));     
end

% set vertical lines at current location and window edges
handles.curLineSr = line([handles.xR(1) handles.xR(1)],[1 nCells]); set(handles.curLineSr,'Color',[0.5 0.5 0.5]);
handles.fwdLineSr = line([handles.xR(1)+0.1 handles.xR(1)+0.1],[1 nCells]); set(handles.fwdLineSr,'Color',[0.3 0.3 0.3]);
handles.revLineSr = line([handles.xR(1)-0.1 handles.xR(1)-0.1],[1 nCells]); set(handles.revLineSr,'Color',[0.3 0.3 0.3]);


% hc panel (left)
axes(handles.BottomLeftAxes); hold on; axis off;
nCells = length(handles.field_loc_L);
for iF = 1:nCells
   
    spks = Data(handles.pf.sd_orig.S{handles.field_cell_L(iF)});
    plot(spks,handles.field_rank_L(iF).*ones(size(spks)),'.','Color',handles.field_col_L(iF,:));
        
end

% set vertical lines at current location and window edges
handles.curLineL = line([handles.xR(1) handles.xR(1)],[1 nCells]); set(handles.curLineL,'Color',[0.5 0.5 0.5]);
handles.fwdLineL = line([handles.xR(1)+0.1 handles.xR(1)+0.1],[1 nCells]); set(handles.fwdLineL,'Color',[0.3 0.3 0.3]);
handles.revLineL = line([handles.xR(1)-0.1 handles.xR(1)-0.1],[1 nCells]); set(handles.revLineL,'Color',[0.3 0.3 0.3]);

% find where landmarks fall in field distribution for horizontal lines
temp = sort(handles.field_loc_L,'ascend'); temp(temp < handles.pf.l.Markers.StartOfMaze) = Inf;
[val,ind] = min(temp-handles.pf.l.Markers.StartOfMaze);
h = line([handles.xR(1) handles.xR(end)], [ind-0.5 ind-0.5]); set(h,'Color',[0.5 0.5 0.5]);

temp = sort(handles.field_loc_L,'ascend'); temp(temp < handles.pf.l.Markers.ct) = Inf;
[val,ind] = min(temp-handles.pf.l.Markers.ct);
h = line([handles.xR(1) handles.xR(end)], [ind-0.5 ind-0.5]); set(h,'Color',[0.5 0.5 0.5]);

temp = sort(handles.field_loc_L,'ascend'); temp(temp < handles.pf.l.Markers.cp) = Inf;
[val,ind] = min(temp-handles.pf.l.Markers.cp);
h = line([handles.xR(1) handles.xR(end)], [ind-0.5 ind-0.5]); set(h,'Color',[0.5 0.5 0.5]);

temp = sort(handles.field_loc_L,'ascend'); temp(temp < handles.pf.l.Markers.f1) = Inf;
[val,ind] = min(temp-handles.pf.l.Markers.f1);
h = line([handles.xR(1) handles.xR(end)], [ind-0.5 ind-0.5]); set(h,'Color',[0.5 0.5 0.5]);

temp = sort(handles.field_loc_L,'ascend'); temp(temp < handles.pf.l.Markers.f2) = Inf;
[val,ind] = min(temp-handles.pf.l.Markers.f2);
h = line([handles.xR(1) handles.xR(end)], [ind-0.5 ind-0.5]); set(h,'Color',[0.5 0.5 0.5]);

% hc panel (right)
axes(handles.BottomRightAxes); hold on; axis off;
for iF = 1:length(handles.field_loc_R)
   
    spks = Data(handles.pf.sd_orig.S{handles.field_cell_R(iF)});
    plot(spks,handles.field_rank_R(iF).*ones(size(spks)),'.','Color',handles.field_col_R(iF,:));
        
end

handles.curLineR = line([handles.xR(1) handles.xR(1)],[1 nCells]); set(handles.curLineR,'Color',[0.5 0.5 0.5]);
handles.fwdLineR = line([handles.xR(1)+0.1 handles.xR(1)+0.1],[1 nCells]); set(handles.fwdLineR,'Color',[0.3 0.3 0.3]);
handles.revLineR = line([handles.xR(1)-0.1 handles.xR(1)-0.1],[1 nCells]); set(handles.revLineR,'Color',[0.3 0.3 0.3]);

% find where landmarks fall in field distribution
temp = sort(handles.field_loc_R,'ascend'); temp(temp < handles.pf.r.Markers.StartOfMaze) = Inf;
[val,ind] = min(temp-handles.pf.r.Markers.StartOfMaze);
h = line([handles.xR(1) handles.xR(end)], [ind-0.5 ind-0.5]); set(h,'Color',[0.5 0.5 0.5]);

temp = sort(handles.field_loc_R,'ascend'); temp(temp < handles.pf.r.Markers.ct) = Inf;
[val,ind] = min(temp-handles.pf.r.Markers.ct);
h = line([handles.xR(1) handles.xR(end)], [ind-0.5 ind-0.5]); set(h,'Color',[0.5 0.5 0.5]);

temp = sort(handles.field_loc_R,'ascend'); temp(temp < handles.pf.r.Markers.cp) = Inf;
[val,ind] = min(temp-handles.pf.r.Markers.cp);
h = line([handles.xR(1) handles.xR(end)], [ind-0.5 ind-0.5]); set(h,'Color',[0.5 0.5 0.5]);

temp = sort(handles.field_loc_R,'ascend'); temp(temp < handles.pf.r.Markers.f1) = Inf;
[val,ind] = min(temp-handles.pf.r.Markers.f1);
h = line([handles.xR(1) handles.xR(end)], [ind-0.5 ind-0.5]); set(h,'Color',[0.5 0.5 0.5]);

temp = sort(handles.field_loc_R,'ascend'); temp(temp < handles.pf.r.Markers.f2) = Inf;
[val,ind] = min(temp-handles.pf.r.Markers.f2);
h = line([handles.xR(1) handles.xR(end)], [ind-0.5 ind-0.5]); set(h,'Color',[0.5 0.5 0.5]);

% Update handles structure
guidata(hObject, handles);

function handles = centerFieldCells(hObject,handles,t)

t = TimeToInd(t,handles.xR);

axes(handles.LeftStr);
axis([handles.xR(t)-0.5 handles.xR(t)+0.5 1 length(handles.pf.sd_str.S)]);
set(handles.curLineSl,'XData',[handles.xR(t) handles.xR(t)]);
set(handles.fwdLineSl,'XData',[handles.xR(t)+0.1 handles.xR(t)+0.1]);
set(handles.revLineSl,'XData',[handles.xR(t)-0.1 handles.xR(t)-0.1]);

axes(handles.RightStr);
axis([handles.xR(t)-0.5 handles.xR(t)+0.5 1 length(handles.pf.sd_str.S)]);
set(handles.curLineSr,'XData',[handles.xR(t) handles.xR(t)]);
set(handles.fwdLineSr,'XData',[handles.xR(t)+0.1 handles.xR(t)+0.1]);
set(handles.revLineSr,'XData',[handles.xR(t)-0.1 handles.xR(t)-0.1]);

axes(handles.BottomLeftAxes);
axis([handles.xR(t)-0.5 handles.xR(t)+0.5 1 length(handles.field_cell_L)]);
set(handles.curLineL,'XData',[handles.xR(t) handles.xR(t)]);
set(handles.fwdLineL,'XData',[handles.xR(t)+0.1 handles.xR(t)+0.1]);
set(handles.revLineL,'XData',[handles.xR(t)-0.1 handles.xR(t)-0.1]);

axes(handles.BottomRightAxes);
axis([handles.xR(t)-0.5 handles.xR(t)+0.5 1 length(handles.field_cell_R)]);
set(handles.curLineR,'XData',[handles.xR(t) handles.xR(t)]);
set(handles.fwdLineR,'XData',[handles.xR(t)+0.1 handles.xR(t)+0.1]);
set(handles.revLineR,'XData',[handles.xR(t)-0.1 handles.xR(t)-0.1]);

% Update handles structure
guidata(hObject, handles);



function CPedit_Callback(hObject, eventdata, handles)
% hObject    handle to CPedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CPedit as text
%        str2double(get(hObject,'String')) returns contents of CPedit as a double


% --- Executes during object creation, after setting all properties.
function CPedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CPedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CPgo.
function CPgo_Callback(hObject, eventdata, handles)
% hObject    handle to CPgo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cpno = str2double(get(handles.CPedit,'String'));
t = handles.pf.cp_entries(cpno);
%[val,t] = min(abs(handles.xR-t));

set(handles.CurrentEdit,'String',num2str(t));

% update rat location
handles = updateRat(hObject,handles,t);
handles = plotMazeCells(hObject,handles,t);
handles = centerFieldCells(hObject,handles,t);
handles = plotRecon(hObject,handles,t);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in CPrev.
function CPrev_Callback(hObject, eventdata, handles)
% hObject    handle to CPrev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cpno = str2double(get(handles.CPedit,'String'));
cpno = max(cpno - 1,1);
t = handles.pf.cp_entries(cpno);
%[val,t] = min(abs(handles.xR-t));

set(handles.CPedit,'String',num2str(cpno));
set(handles.CurrentEdit,'String',num2str(t));

% update rat location
t = TimeToInd(t,handles.xR);
handles = updateRat(hObject,handles,t);
handles = plotMazeCells(hObject,handles,t);
handles = centerFieldCells(hObject,handles,t);
handles = plotRecon(hObject,handles,t);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in CPfwd.
function CPfwd_Callback(hObject, eventdata, handles)
% hObject    handle to CPfwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cpno = str2double(get(handles.CPedit,'String'));
cpno = min(cpno + 1,str2double(get(handles.CPmax,'String')));
t = handles.pf.cp_entries(cpno);
%[val,t] = min(abs(handles.xR-t));

set(handles.CPedit,'String',num2str(cpno));
set(handles.CurrentEdit,'String',num2str(t));

% update rat location
handles = updateRat(hObject,handles,t);
handles = plotMazeCells(hObject,handles,t);
handles = centerFieldCells(hObject,handles,t);
handles = plotRecon(hObject,handles,t);

% Update handles structure
guidata(hObject, handles);

function ind = TimeToInd(t,tvec)
% return index of tvec element closest to input time t

[val,ind] = min(abs(tvec-t));

function handles = plotRecon(hObject,handles,t)

t = TimeToInd(t,handles.xR);
pass = handles.cpp(t);
if ~isnan(pass) % in cp pass
    
    if isfield(handles,'cp')
    
    [val,ind] = min(abs(handles.cp.tvec{pass}-handles.xR(t)));
    axes(handles.TopCenterAxes); cla;
    temp = nan(63,47); temp(handles.cp.goodOccInd) = handles.cp.p5{pass}(ind,:);
    h = imagesc(temp'); axis xy; hold on;
    %plot(handles.cp.yBinned{pass}(ind),handles.cp.xBinned{pass}(ind),'xw','MarkerSize',30);
    plot(handles.cp.xBinned{pass}(ind),handles.cp.yBinned{pass}(ind),'xw','MarkerSize',30);
    %set(gca,'YDir','reverse'); set(gca,'XDir','reverse');
    axis off;
    
    axes(handles.BottomCenterAxes); cla;
    temp = nan(63,47); temp(handles.cp.goodOccInd) = handles.cp.p0{pass}(ind,:);
    h = imagesc(temp'); axis xy; hold on;
    %plot(handles.cp.yBinned{pass}(ind),handles.cp.xBinned{pass}(ind),'xw','MarkerSize',30);
    plot(handles.cp.xBinned{pass}(ind),handles.cp.yBinned{pass}(ind),'xw','MarkerSize',30);
    %set(gca,'YDir','reverse'); set(gca,'XDir','reverse');
    axis off;
    
    end
    
end

% Update handles structure
guidata(hObject, handles);



function swedit_Callback(hObject, eventdata, handles)
% hObject    handle to swedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of swedit as text
%        str2double(get(hObject,'String')) returns contents of swedit as a double


% --- Executes during object creation, after setting all properties.
function swedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to swedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SWgo.
function SWgo_Callback(hObject, eventdata, handles)
% hObject    handle to CPgo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

swno = str2double(get(handles.swedit,'String'));
t = handles.pf.swt(swno);
%[val,t] = min(abs(handles.xR-t));

set(handles.CurrentEdit,'String',num2str(t));

% update rat location
handles = updateRat(hObject,handles,t);
handles = plotMazeCells(hObject,handles,t);
handles = centerFieldCells(hObject,handles,t);
handles = plotRecon(hObject,handles,t);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in SWrev.
function SWrev_Callback(hObject, eventdata, handles)
% hObject    handle to SWrev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SWfwd.
function SWfwd_Callback(hObject, eventdata, handles)
% hObject    handle to CPfwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

swno = str2double(get(handles.swedit,'String'));
swno = min(swno + 1,str2double(get(handles.SWmax,'String')));
t = handles.pf.swt(swno);
%[val,t] = min(abs(handles.xR-t));

set(handles.swedit,'String',num2str(swno));
set(handles.CurrentEdit,'String',num2str(t));

% update rat location
handles = updateRat(hObject,handles,t);
handles = plotMazeCells(hObject,handles,t);
handles = centerFieldCells(hObject,handles,t);
handles = plotRecon(hObject,handles,t);

% Update handles structure
guidata(hObject, handles);
