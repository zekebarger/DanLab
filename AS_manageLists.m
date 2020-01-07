% zeke barger 120619
% TODO disable addbelow... it isn't getting updated :/

function varargout = AS_manageLists(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AS_manageLists_OpeningFcn, ...
    'gui_OutputFcn',  @AS_manageLists_OutputFcn, ...
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


% --- Executes just before AS_manageLists is made visible.
function AS_manageLists_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AS_manageLists (see VARARGIN)

% Choose default command line output for AS_manageLists
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Store data in an invisible text box called C
% it's not great, but beats using a global variable
setappdata(handles.C,'line',1);
setappdata(handles.C,'text',{});
setappdata(handles.C,'dirlistIsCurrent',0); % is the dirlist being displayed equal to
%                                             the saved file
setappdata(handles.C,'current1',{}); % current state of the dirlist
setappdata(handles.C,'current2',{}); % current state of the mouselist
setappdata(handles.C,'history1',{}); % past states of the dirlist
setappdata(handles.C,'history2',{}); % past states of the mouselist
setappdata(handles.C,'future1',{}); % "future" states of the dirlist
setappdata(handles.C,'future2',{}); % "future" states of the mouselist 
setappdata(handles.C,'lastInput1',''); % last path/file chosen
setappdata(handles.C,'lastOutput1',''); % last path/file chosen
setappdata(handles.C,'lastInput2',''); 
setappdata(handles.C,'lastOutput2',''); 
setappdata(handles.C,'listboxes',[handles.listbox1,handles.listbox2]);

set(handles.text2,'HorizontalAlignment','left')
% clear boxes
set(handles.text2,'String',{});
set(handles.listbox1,'string',{});
set(handles.listbox1,'Value',1);
set(handles.listbox2,'string',{});
set(handles.listbox2,'Value',1);
toggle_dirlistCurrent(handles,0);

% --- Outputs from this function are returned to the command line.
function varargout = AS_manageLists_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function disptext(handles, string)
txt = getappdata(handles.C,'text');
txt{getappdata(handles.C,'line')} = string;
setappdata(handles.C,'text',txt);
setappdata(handles.C,'line',getappdata(handles.C,'line')+1);

if getappdata(handles.C,'line') == 15
    setappdata(handles.C,'line', 14)
    txt = getappdata(handles.C,'text');
    setappdata(handles.C,'text',txt(2:14));
end
set(handles.text2,'String',getappdata(handles.C,'text'));

function backup(handles, p)
ps = num2str(p);
history = getappdata(handles.C,['history',ps]);
if isempty(history)
    history = {{}};
else
    history{end+1} = getappdata(handles.C,['current',ps]);
end
setappdata(handles.C,['history',ps],history);
setappdata(handles.C,['future',ps],{});
if ps==1
    toggle_dirlistCurrent(handles,0);
end

function alphaFcn(handles,p)
ps = num2str(p);
backup(handles, p)
setappdata(handles.C,['current',ps],sort(getappdata(handles.C,['current',ps])))
listboxes = getappdata(handles.C,'listboxes');
set(listboxes(p),'Value',1);
set(listboxes(p),'string',getappdata(handles.C,['current',ps]));

function alphabtn1_Callback(hObject, eventdata, handles)
alphaFcn(handles,1)

function alphabtn2_Callback(hObject, eventdata, handles)
alphaFcn(handles,2)

function savebtn1_Callback(hObject, eventdata, handles)
[FileName,PathName,~] = uiputfile('*','save dirlist',...
    getStringUntilLastSlash(getappdata(handles.C,'lastOutput1')));
if ischar(FileName) % user didn't cancel
    if ~contains(FileName,'.mat')
        FileName = [FileName,'.mat'];
    end
    setappdata(handles.C,'lastOutput1',[PathName,FileName]);
    dirlist = getappdata(handles.C,'current1');
    save(cat(2,PathName,FileName),'dirlist');
    disptext(handles,cat(2,'saved ',FileName));
    toggle_dirlistCurrent(handles,1);
end

function savebtn2_Callback(hObject, eventdata, handles)
[FileName,PathName,~] = uiputfile('*','save mouselist',...
    getappdata(handles.C,'lastOutput2'));
if ischar(FileName) % user didn't cancel
    if ~contains(FileName,'.mat')
        FileName = [FileName,'.mat'];
    end
    setappdata(handles.C,'lastOutput2',getStringUntilLastSlash(PathName));
    mouselist = getappdata(handles.C,'current2');
    save(cat(2,PathName,FileName),'mouselist');
    disptext(handles,cat(2,'saved ',FileName));
end

function addFcn(handles,p)
ps = num2str(p);
if p==1
    PathName = uigetdir(getappdata(handles.C,'lastInput1'),'Choose recording folder');
else
    [FileName,PathName,~] = uigetfile('*.mat','Choose mouselist file',...
        getappdata(handles.C,'lastInput2'));
    PathName = [PathName,FileName];
end
% check if something was selected
if ~ischar(PathName)
    return
end
% store the path for convenience
setappdata(handles.C,['lastInput',ps],getStringUntilLastSlash(PathName));

% check if the entry is already present
listX = getappdata(handles.C,['current',ps]);
if any(strcmp(PathName,listX))
    disptext(handles,'entry is already present');
    return
end

if p==2 % check if the file actually contains a dirlist
    if ~containsVariable(PathName,'dirlist')
        disptext(handles,'file does not contain a dirlist');
        return
    end
end

backup(handles,p);
setappdata(handles.C,['current',ps],cat(2,getappdata(handles.C,['current',ps]),PathName));
listboxes = getappdata(handles.C,'listboxes');
set(listboxes(p),'string',getappdata(handles.C,['current',ps]));
disptext(handles,'added entry');

function addbtn1_Callback(hObject, eventdata, handles)
addFcn(handles,1)

function addbtn2_Callback(hObject, eventdata, handles)
addFcn(handles,2)

function removeFcn(handles,p)
ps = num2str(p);
listX = getappdata(handles.C,['current',ps]);
if isempty(listX)
    return
end
backup(handles, p);
listboxes = getappdata(handles.C,'listboxes');
index_selected = get(listboxes(p),'Value');
disptext(handles,cat(2,'removed ',listX{index_selected}));
listX(index_selected) = [];
setappdata(handles.C,['current',ps], listX)
set(listboxes(p),'string',listX);
set(listboxes(p),'value',max([min([index_selected,length(listX)]),1]));

function removebtn1_Callback(hObject, eventdata, handles)
removeFcn(handles,1)

function removebtn2_Callback(hObject, eventdata, handles)
removeFcn(handles,2)

function subtractFcn(handles,p)
ps = num2str(p);
[FileName,PathName,~] = uigetfile('*.mat','Choose file to subtract',...
    getappdata(handles.C,['lastInput',ps]));
% check if something was selected
if ~ischar(FileName)
    disptext(handles,'no file selected');
    return
end
% store location for convenience
setappdata(handles.C,['lastInput',ps],getStringUntilLastSlash(PathName));

% check if file contains the correct variable
if p==1
	variableName = 'dirlist';
else
    variableName = 'mouselist';
end
if ~containsVariable([PathName,FileName], variableName)
    disptext(handles,['file does not contain a ', variableName]);
    return
end

temp = load(cat(2,PathName,FileName));
if p==1
    list2sub=temp.dirlist;
else
    list2sub=temp.mouselist;
end
disptext(handles,cat(2,'loaded ',FileName));
if isempty(list2sub)
    disptext(handles,'but the list to be subtracted is empty');
    return
end

% check for entries not in our current list
listX = getappdata(handles.C,['current',ps]);
notInList = zeros(1,length(list2sub));
for i = 1:length(list2sub)
    notInList(i) = ~any(strcmp(list2sub{i},listX));
end
% remove those from the list to be subtracted
list2sub(logical(notInList)) = [];
if isempty(list2sub)
    return
end

% if there's something to subtract
backup(handles,p);
listX = getappdata(handles.C,['current',ps]);
indicesToRemove = zeros(1,length(listX));
for i = 1:length(listX)
    indicesToRemove(i) = any(strcmp(listX{i},list2sub));
end
listX(logical(indicesToRemove))=[];
setappdata(handles.C,['current',ps],listX)

listboxes = getappdata(handles.C,'listboxes');
set(listboxes(p),'string',getappdata(handles.C,['current',ps]));
set(listboxes(p),'value',length(getappdata(handles.C,['current',ps])));
disptext(handles,['Removed ', num2str(length(list2sub)),' entries',...
    ', ignored ',num2str(sum(notInList)),' non-overlapping entries']);

function subtractbtn1_Callback(hObject, eventdata, handles)
subtractFcn(handles,1)

function subtractbtn2_Callback(hObject, eventdata, handles)
subtractFcn(handles,2)

function loadFcn(handles, p)
ps = num2str(p);
if p==1
    variableName = 'dirlist';
else
    variableName = 'mouselist';
end
[FileName,PathName,~] = uigetfile('*.mat',['Choose file containing a ',variableName],...
    getappdata(handles.C,['lastInput',ps]));
% check that something was selected
if ~ischar(FileName)
    disptext(handles,'no file selected');
    return
end
% store location for convenience
setappdata(handles.C,['lastInput',ps],getStringUntilLastSlash(PathName));

% check if file contains the correct type of list
if ~containsVariable([PathName,FileName],variableName)
    disptext(handles,['file does not contain a ',variableName]);
    return
end

backup(handles, p);
temp = load([PathName,FileName]);
if p==1
    setappdata(handles.C,'current1',temp.dirlist);
    toggle_dirlistCurrent(handles,1);
    setappdata(handles.C,'lastOutput1',[PathName,FileName]);
else
    setappdata(handles.C,'current2',temp.mouselist);
end
listboxes = getappdata(handles.C,'listboxes');
set(listboxes(p),'string',getappdata(handles.C,['current',ps]));
set(listboxes(p),'value',1);
disptext(handles,cat(2,'loaded ',FileName));

function loadbtn1_Callback(hObject, eventdata, handles)
loadFcn(handles, 1)

function loadbtn2_Callback(hObject, eventdata, handles)
loadFcn(handles, 2)

function clearFcn(handles, p)
ps = num2str(p);
if isempty(getappdata(handles.C,['current',ps]))
    return
end
backup(handles, p);
setappdata(handles.C,['current',ps],{})
listboxes = getappdata(handles.C,'listboxes');
set(listboxes(p),'Value',1);
set(listboxes(p),'string',getappdata(handles.C,['current',ps]));

function clearbtn1_Callback(hObject, eventdata, handles)
clearFcn(handles, 1)

function clearbtn2_Callback(hObject, eventdata, handles)
clearFcn(handles, 2)

function appendFcn(handles,p)
ps = num2str(p);
if p==1
    variableName = 'dirlist';
else
    variableName = 'mouselist';
end
[FileName,PathName,~] = uigetfile('*.mat',['Choose file containing a ',variableName],...
    getappdata(handles.C,['lastInput',ps]));
% check if something was selected
if ~ischar(FileName)
    disptext(handles,'no file selected');
    return
end
% store location for convenience
setappdata(handles.C,['lastInput',ps],getStringUntilLastSlash(PathName));

% check if file contains the correct type of list
if ~containsVariable([PathName,FileName],variableName)
    disptext(handles,['file does not contain a ', variableName]);
    return
end

temp = load(cat(2,PathName,FileName));
if p==1
    list2add=temp.dirlist;
else
    list2add=temp.mouselist;
end
disptext(handles,cat(2,'loaded ',FileName));

% check for entries that are already present
listX = getappdata(handles.C,['current',ps]);
alreadyInList = zeros(1,length(list2add));
for i = 1:length(list2add)
    alreadyInList(i) = any(strcmp(list2add{i},listX));
end
% remove those from the list to be added
list2add(logical(alreadyInList)) = [];
if isempty(list2add)
    return
end

% if there's something to add
backup(handles, p);
setappdata(handles.C,['current',ps],[getappdata(handles.C,['current',ps]),list2add])
listboxes = getappdata(handles.C,'listboxes');
set(listboxes(p),'value',1);
set(listboxes(p),'string',getappdata(handles.C,['current',ps]));
set(listboxes(p),'value',length(getappdata(handles.C,['current',ps])));
disptext(handles,['Added ', num2str(length(list2add)),' entries(s)',...
    ', ignored ',num2str(sum(alreadyInList)),' duplicate(s)']);


function appendbtn1_Callback(hObject, eventdata, handles)
appendFcn(handles,1)

function appendbtn2_Callback(hObject, eventdata, handles)
appendFcn(handles,2)

function undoFcn(handles,p)
ps = num2str(p);
historyVar = ['history',ps];
currentVar = ['current',ps];
futureVar = ['future',ps];
% if there's no history to return to, don't do anything
history = getappdata(handles.C,historyVar);
if isempty(history)
    return
end

% store current state in the beginning of the future list
future = getappdata(handles.C,futureVar);
if isempty(future)
    setappdata(handles.C,futureVar,{getappdata(handles.C,currentVar)})
else    
    setappdata(handles.C,futureVar,[{getappdata(handles.C,currentVar)},future]);
end

% take the most recent item in history and set it to the present
setappdata(handles.C,currentVar,history{end});

% remove the end of history
history(end) = [];
setappdata(handles.C,historyVar, history);

% update display
listboxes = getappdata(handles.C,'listboxes');
set(listboxes(p),'string',getappdata(handles.C,currentVar));
set(listboxes(p),'value',length(getappdata(handles.C,currentVar)));
disptext(handles,'undo');

if p==1
    toggle_dirlistCurrent(handles,0);
end

function undobtn1_Callback(hObject, eventdata, handles)
undoFcn(handles,1)

function undobtn2_Callback(hObject, eventdata, handles)
undoFcn(handles,2)

function redoFcn(handles,p)
ps = num2str(p);
historyVar = ['history',ps];
currentVar = ['current',ps];
futureVar = ['future',ps];
% if there's no future to advance to, don't do anything
future = getappdata(handles.C,futureVar);
if isempty(future)
    return
end

% store current state in the end of history
history = getappdata(handles.C,historyVar);
if isempty(history)
    setappdata(handles.C,historyVar,{getappdata(handles.C,currentVar)})
else    
    setappdata(handles.C,historyVar,[history,{getappdata(handles.C,currentVar)}]);
end

% take the most recent item in future and set it to the present
setappdata(handles.C,currentVar,future{1});

% remove the beginning of the future
future(1) = [];
setappdata(handles.C,futureVar, future);

% update display
listboxes = getappdata(handles.C,'listboxes');
set(listboxes(p),'string',getappdata(handles.C,currentVar));
set(listboxes(p),'value',length(getappdata(handles.C,currentVar)));
disptext(handles,'redo');

if p==1
    toggle_dirlistCurrent(handles,0);
end

function redobtn1_Callback(hObject, eventdata, handles)
redoFcn(handles,1)

function redobtn2_Callback(hObject, eventdata, handles)
redoFcn(handles,2)

function driveFcn(handles,p)
ps = num2str(p);
listX = getappdata(handles.C,['current',ps]);
if isempty(listX)
    return
end

% open dialog box to get new drive name
idx = strfind(listX{1},'\');
oldName = listX{1}(1:(idx(1)-1));
name = inputdlg('Enter new hard drive name for all entries in the list (including the colon):',...
    'Drive name',[1 40],{oldName});
% if something was entered
if ~isempty(name)
    name = name{1};
    backup(handles,2);
    % remove any slashes from beginning or end
    name = strip(name,'\');
    name = strip(name,'/');
    % replace top level folder in all entries
    for i = 1:length(listX)
        % get the string
        s = listX{i};
        % remove up to first slash
        idx = strfind(s,'\');
        listX{i} = [name,s(idx(1):end)];
    end
    setappdata(handles.C,['current',ps],listX);
    % update display
    listboxes = getappdata(handles.C,'listboxes');
    set(listboxes(p),'string',getappdata(handles.C,['current',ps]));
    disptext(handles,'replaced drive name');
end

function drivebtn1_Callback(hObject, eventdata, handles)
driveFcn(handles,1)

function drivebtn2_Callback(hObject, eventdata, handles)
driveFcn(handles,2)

% check if file on some path contains a variable with some name 
function [present] = containsVariable(path,varName)
fileContents = whos('-file',path); % get information about file contents
varInfo = []; % just keep info about the named variable (if it exists)
for i = 1:length(fileContents)
    if strcmp(varName,fileContents(i).name)
        varInfo = fileContents(i);
        break
    end
end
present = ~isempty(varInfo);

function mergebtn_Callback(hObject, eventdata, handles)
if length(getappdata(handles.C, 'current2')) < 2
    warndlg('At least two lists must be present');
    return
end

[FileName,PathName,~] = uiputfile('*','save merged mouselists as dirlist',...
    getappdata(handles.C,'lastOutput2'));
if ischar(FileName) % user didn't cancel
    % store location for convenience
    setappdata(handles.C,'lastOutput2',getStringUntilLastSlash(PathName));
    dirlist = {};
    mouselist = getappdata(handles.C,'current2');
    for i = 1:length(mouselist)
        d = load(mouselist{i});
        dirlist = [dirlist, d.dirlist];
    end
    save(cat(2,PathName,FileName),'dirlist');
    disptext(handles,cat(2,'saved ',FileName));
end

% --- Executes on button press in addbelowbtn.
function addbelowbtn_Callback(hObject, eventdata, handles)
% check if entry is already present
mouselist = getappdata(handles.C,'current2');
for j = 1:(numel(mouselist)) % go through old list
    if strcmp(mouselist{j},getappdata(handles.C,'lastOutput1')) % same entry
        disptext(handles,'entry is already present');
        return
    end
end

% check if file actually contains a dirlist
if ~containsVariable(getappdata(handles.C,'lastOutput1'),'dirlist')
    disptext(handles,'file does not contain a dirlist');
    return
end

backup(handles,2);
setappdata(handles.C,'current2',...
    cat(2,getappdata(handles.C,'current2'),getappdata(handles.C,'lastOutput1')));
set(handles.listbox2,'string',getappdata(handles.C,'current2'));
set(handles.listbox2,'value',...
    length(getappdata(handles.C,'current2')));
disptext(handles,'added entry');

function toggle_dirlistCurrent(handles,current)
setappdata(handles.C,'dirlistIsCurrent',current);
if current
    set(handles.addbelowbtn,'Enable','on');
else
    set(handles.addbelowbtn,'Enable','off');
end

function [s] = getStringUntilLastSlash(s)
if isempty(s)
    return
end
idx = strfind(s,'\');
s = s(1 : (idx(end)-1));

% --- Executes on key release with focus on figure1 or any of its controls.
function figure1_WindowKeyReleaseFcn(hObject, eventdata, handles)

function listbox1_Callback(hObject, eventdata, handles)

function listbox1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox2_Callback(hObject, eventdata, handles)

function listbox2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
