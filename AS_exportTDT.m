function varargout = AS_exportTDT(varargin)
% zeke barger 081120
% Exports EEG/EMG recording from TDT system into a format that can be
% interpreted by AccuSleep. The onsets and offsets of laser stimuli, if
% present are saved (in seconds).

% TODO: don't overwrite folders when making them

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AS_exportTDT_OpeningFcn, ...
    'gui_OutputFcn',  @AS_exportTDT_OutputFcn, ...
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


% --- Executes just before AS_exportTDT is made visible.
function AS_exportTDT_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for AS_exportTDT
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% Store data in an invisible text box called C
% it's not great, but beats using a global variable
setappdata(handles.C,'blockPath',[]);
setappdata(handles.C,'outputLocation1',[]);
setappdata(handles.C,'outputLocation2',[]);
setappdata(handles.C,'outputName1',[]);
setappdata(handles.C,'outputName2',[]);
setappdata(handles.C,'SR1',512);
setappdata(handles.C,'SR2',512);
setappdata(handles.C,'EEGch1',[]);
setappdata(handles.C,'EMGch1',[]);
setappdata(handles.C,'EEGch2',[]);
setappdata(handles.C,'EMGch2',[]);

% console text
setappdata(handles.messagebox,'text',{});
setappdata(handles.messagebox,'line',1);

% clear all fields
set(handles.blockbox,'String','');
set(handles.eegbox,'String','');
set(handles.emgbox,'String','');
set(handles.srbox,'String','512');
set(handles.locationbox,'String','');
set(handles.namebox,'String','');
set(handles.eegbox2,'String','');
set(handles.emgbox2,'String','');
set(handles.srbox2,'String','512');
set(handles.locationbox2,'String','');
set(handles.namebox2,'String','');
set(handles.messagebox,'String',{});

function varargout = AS_exportTDT_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function blockbtn_Callback(~, ~, handles)
% choose the block
blockPath = uigetdir(get(handles.blockbox,'String'), 'Choose block to process');

% make sure something was selected
if length(blockPath) < 2
    return
end

% update the display and store the path
setappdata(handles.C,'blockPath',blockPath);
set(handles.blockbox,'String',blockPath);

function blockbox_Callback(hObject, ~, handles)
% update the stored path
setappdata(handles.C,'blockPath',get(hObject,'String'));

function blockbox_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function previewbtn_Callback(~, ~, handles)
% determine if path is valid
blockPath = getappdata(handles.C,'blockPath');
if ~exist(blockPath,'dir')
    
    ('No block selected');
    return
end

% load the data, if it's not already loaded
setListLock(handles,1);

try
    disptext(handles, 'Loading, please wait...')
    drawnow;
    % collect basic info about the block
    basic_data = TDT2mat(blockPath,'VERBOSE',0,'T1',0,'T2',1);
    durationInSeconds = seconds(duration(basic_data.info.duration));
    % decide how much data to load (20 mins tops)
    T2 = min([20*60, durationInSeconds]);
    % actually load the data
    blockData = TDT2mat(blockPath,'VERBOSE',0,'T1',0,'T2',T2);
catch
    errordlg(['Could not load data from the block. Either the block ',...
        'is not properly formatted, or you need to install OpenDeveloper.']);
    setListLock(handles,0);
    return
end

disptext(handles, 'Finished loading.')
drawnow;

% get information about the data
numChannels = size(blockData.streams.LFPs.data,1);
numSamples = size(blockData.streams.LFPs.data,2);
SR = blockData.streams.LFPs.fs;
% calculate min, max, text position for each channel
stats = zeros(numChannels, 3);
for i = 1:numChannels
    stats(i,1) = min(blockData.streams.LFPs.data(i,:));
    stats(i,2) = max(blockData.streams.LFPs.data(i,:));
    stats(i,3) = (stats(i,2) + stats(i,1))/2;
end
% if some channels are all 0, make sure to space them properly
ranges = stats(:,2)-stats(:,1);
if ~any(ranges) % case where all channels are just 0
    minHeight = 1;
else
    minHeight = min(ranges(ranges>0));
end
for i = 1:numChannels
    if ranges(i)==0
        stats(i,1) = 0;
        stats(i,2) = minHeight;
        stats(i,3) = 0;
    end
end
warning('off','MATLAB:colon:nonIntegerIndex');

labels = cell(1,numChannels);
figure
plot(blockData.streams.LFPs.data(1,1:(SR/200):end));
labels{1} = 'Channel 1';
for i = 2:numChannels
    hold on, plot(blockData.streams.LFPs.data(i,1:(SR/200):end)-...
        (stats(i,2)-stats(i-1,1)));
    stats(i,:) = stats(i,:) - (stats(i,2)-stats(i-1,1));
    labels{i} = ['Channel ',num2str(i)];
end
axis tight
title('Preview (~20 mins)')
set(gca,'YTick',flipud(stats(:,3)),'YTickLabel',fliplr(labels),'XTick',[])
warning('on','MATLAB:colon:nonIntegerIndex');
setListLock(handles,0);

function locationbtn_Callback(~, ~, handles)
% choose the location for the output folder
location = uigetdir(get(handles.locationbox,'String'), 'Choose location for the output folder');

% make sure something was selected
if length(location) < 2
    return
end

% update the display and store the path
setappdata(handles.C,'outputLocation1',location);
set(handles.locationbox,'String',location);

function locationbox_Callback(hObject, ~, handles)
% update the stored location
setappdata(handles.C,'outputLocation1',get(hObject,'String'));

function locationbox_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function namebox_Callback(hObject, ~, handles)
% update the stored name
s = get(hObject,'String');
s = strip(s,'\');
s = strip(s,'/');
set(handles.namebox,'String',s)
setappdata(handles.C,'outputName1',s);

function namebox_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function srbox_Callback(hObject, ~, handles)
% update the stored SR
setappdata(handles.C,'SR1',str2num(get(hObject,'String')));

function srbox_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function eegbox_Callback(hObject, ~, handles)
% update the stored EEG channel
setappdata(handles.C,'EEGch1',str2num(get(hObject,'String')));

function eegbox_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function emgbox_Callback(hObject, ~, handles)
% update the stored EMG channel
setappdata(handles.C,'EMGch1',str2num(get(hObject,'String')));

function emgbox_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function runbtn_Callback(~, ~, handles)
general_run(handles, 1);

% add a line of text to the message box
function disptext(handles, string)
t = getappdata(handles.messagebox,'text');
ln = getappdata(handles.messagebox,'line');

t{ln}=string;
ln = ln+1;
if ln == 9
    ln = 8;
    t = t(2:8);
end
set(handles.messagebox,'String',t)
setappdata(handles.messagebox,'text',t);
setappdata(handles.messagebox,'line',ln);

% --- Executes on button press in runbothbtn.
function runbothbtn_Callback(hObject, eventdata, handles)
general_run(handles, 3);

function emgbox2_Callback(hObject, eventdata, handles)
% update the stored EMG channel
setappdata(handles.C,'EMGch2',str2num(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function emgbox2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function eegbox2_Callback(hObject, eventdata, handles)
% update the stored EEG channel
setappdata(handles.C,'EEGch2',str2num(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function eegbox2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runbtn2.
function runbtn2_Callback(hObject, eventdata, handles)
general_run(handles, 2);

function srbox2_Callback(hObject, eventdata, handles)
% update the stored SR
setappdata(handles.C,'SR2',str2num(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function srbox2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function namebox2_Callback(hObject, eventdata, handles)
% update the stored name
s = get(hObject,'String');
s = strip(s,'\');
s = strip(s,'/');
set(handles.namebox2,'String',s)
setappdata(handles.C,'outputName2',s);


% --- Executes during object creation, after setting all properties.
function namebox2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function locationbox2_Callback(hObject, eventdata, handles)
% update the stored location
setappdata(handles.C,'outputLocation2',get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function locationbox2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in locationbtn2.
function locationbtn2_Callback(hObject, eventdata, handles)
% choose the location for the output folder
location = uigetdir(get(handles.locationbox2,'String'), 'Choose location for the output folder');

% make sure something was selected
if length(location) < 2
    return
end

% update the display and store the path
setappdata(handles.C,'outputLocation2',location);
set(handles.locationbox2,'String',location);

% process data from one or both mice
% mouse: 1=mouse1 only, 2 = mouse2 only, 3 =both
function general_run(handles, mouse)
% decide which inputs to use
switch mouse
    case 1
        a = 1;
        b = 1;
    case 2
        a = 2;
        b = 2;
    case 3
        a = 1;
        b = 2;
end

% check that a block is selected
blockPath = getappdata(handles.C,'blockPath');
if ~exist(blockPath,'dir')
    errordlg('No block selected');
    return
end

setListLock(handles,1);

for n = a:b % for each mouse selected
    % check that all inputs are valid
    EEGch = getappdata(handles.C,['EEGch',num2str(n)]);
    if isempty(EEGch)
        errordlg(['No EEG channel selected for mouse ',num2str(n)]);
        setListLock(handles,0);
        return
    end
    
    EMGch = getappdata(handles.C,['EMGch',num2str(n)]);
    if isempty(EMGch)
        errordlg(['No EMG channel selected for mouse ',num2str(n)]);
        setListLock(handles,0);
        return
    end
    
    location = getappdata(handles.C,['outputLocation',num2str(n)]);
    if ~exist(location,'dir')
        errordlg(['No output folder location set for mouse ',num2str(n)]);
        setListLock(handles,0);
        return
    end
    
    outputName = getappdata(handles.C,['outputName',num2str(n)]);
    if isempty(outputName)
        errordlg(['No output folder name set for mouse ',num2str(n)]);
        setListLock(handles,0);
        return
    end
    
    SR = getappdata(handles.C,['SR',num2str(n)]);
    if isempty(SR)
        errordlg(['No sampling rate set for mouse ',num2str(n)]);
        setListLock(handles,0);
        return
    end
    
    % make sure last character of output folder location is \
    if ~strcmp(location(end),filesep)
        location = [location,filesep];
    end
    % opposite for output folder name
    if contains(outputName,{'\','/'})
        errordlg(['Output folder name should not contain slashes. ',...
            'A good folder name might be the date of the experiment']);
        setListLock(handles,0);
        return
    end
    
    % load the data
    try
        disptext(handles, ['Processing mouse ',num2str(n),', please wait...']);
        drawnow;
        % load eeg 
        eegData = TDT2mat(blockPath,'VERBOSE',0,'CHANNEL',EEGch);
        % load emg
        emgData = TDT2mat(blockPath,'VERBOSE',0,'CHANNEL',EMGch);
    catch
        errordlg(['Could not load data from the block. Either the block ',...
            'is not properly formatted, or you need to install OpenDeveloper.']);
        setListLock(handles,0);
        return
    end
    drawnow;
    
    % ok! now we can do things
    warning('off','MATLAB:colon:nonIntegerIndex');
    
    drawnow;
    % make the folder, if necessary
    if exist([location,outputName],'dir') % check that output folder does not already exist
        answer = questdlg([location,outputName,...
            ' already exists. EEG/EMG/laser contents may be overwritten.'], ...
            'Folder already exists', ...
            'Continue','Exit','Exit');
        switch answer
            case 'Exit'
                setListLock(handles,0);
                return
        end
    else
        mkdir([location,outputName]);
    end
    % get EEG, EMG
    originalSR = eegData.streams.LFPs.fs;
    EEG = eegData.streams.LFPs.data(EEGch,1:(originalSR/SR):end);
    EMG = emgData.streams.LFPs.data(EMGch,1:(originalSR/SR):end);
    % get laser stimulus onsets/offsets
    pulse_onsets = eegData.epocs.LasT.onset; % get laser pulse times
    pulse_offsets = eegData.epocs.LasT.offset;
    % the pulse times sometimes contain... infinity... get rid of that
    pulse_onsets(isinf(pulse_onsets)) = [];
    pulse_offsets(isinf(pulse_offsets)) = [];
    % the last onset and offset are usually (always?) artifacts that appear
    % simultaneously
    % so, if we see that the last onset and offset are very close together,
    % they should probably be removed
    if abs(pulse_onsets(end) - pulse_offsets(end)) < 2
        pulse_onsets = pulse_onsets(1:(end-1)); % remove last entry
        pulse_offsets = pulse_offsets(1:(end-1)); % remove last entry
    end
    % get onsets and offsets for pulse trains
    if isempty(pulse_onsets)
        onsets = [];
        offsets = [];
    else
        % check if first offset is before first onset, and remove if so
        while pulse_offsets(1) < pulse_onsets(1)
            pulse_offsets(1) = [];
        end
        % check if last onset is after last offset, and remove if so
        while pulse_onsets(end) > pulse_offsets(end)
            pulse_onsets(end) = [];
        end
        minISI = 30; % max allowable interval between pulses in a train, in seconds
        onsets = [pulse_onsets(1);pulse_onsets(1+find(diff(pulse_onsets) > minISI))];
        offsets = [pulse_offsets(find(diff(pulse_offsets) > minISI)); pulse_offsets(end)];
    end
    
%     % old way of saving laser data 
%     laser = zeros(1,size(blockData.streams.LFPs.data,2)); % preallocate
%     for i = 1:length(onsets)
%         laser(round(originalSR*onsets(i)) : round(originalSR*offsets(i))) = 1;
%     end
%     laser = laser(1:(originalSR/SR):end);

    % new way of saving laser data (more compact)
    if ~isempty(onsets) && ~isempty(offsets)   
        laser = [onsets, offsets];
        save([location,outputName,filesep,'laser.mat'],'laser');
    end
    
    save([location,outputName,filesep,'EEG.mat'],'EEG');
    save([location,outputName,filesep,'EMG.mat'],'EMG');
    
    warning('on','MATLAB:colon:nonIntegerIndex');
    disptext(handles, ['Finished processing data for mouse ',num2str(n)])
    drawnow;
end
setListLock(handles,0);

% lock the list box while the program is busy (1), or unlock it (0)
function setListLock(handles,locked)
if locked
    set(handles.blockbtn,'Enable','off');
    set(handles.blockbox,'Enable','off');
    set(handles.previewbtn,'Enable','off');
    set(handles.eegbox,'Enable','off');
    set(handles.emgbox,'Enable','off');
    set(handles.srbox,'Enable','off');
    set(handles.locationbtn,'Enable','off');
    set(handles.locationbox,'Enable','off');
    set(handles.namebox,'Enable','off');
    set(handles.runbtn,'Enable','off');  
    set(handles.eegbox2,'Enable','off');
    set(handles.emgbox2,'Enable','off');
    set(handles.srbox2,'Enable','off');
    set(handles.locationbtn2,'Enable','off');
    set(handles.locationbox2,'Enable','off');
    set(handles.namebox2,'Enable','off');
    set(handles.runbtn2,'Enable','off');
    set(handles.runbothbtn,'Enable','off');
else
    set(handles.blockbtn,'Enable','on');
    set(handles.blockbox,'Enable','on');
    set(handles.previewbtn,'Enable','on');
    set(handles.eegbox,'Enable','on');
    set(handles.emgbox,'Enable','on');
    set(handles.srbox,'Enable','on');
    set(handles.locationbtn,'Enable','on');
    set(handles.locationbox,'Enable','on');
    set(handles.namebox,'Enable','on');
    set(handles.runbtn,'Enable','on');
    set(handles.eegbox2,'Enable','on');
    set(handles.emgbox2,'Enable','on');
    set(handles.srbox2,'Enable','on');
    set(handles.locationbtn2,'Enable','on');
    set(handles.locationbox2,'Enable','on');
    set(handles.namebox2,'Enable','on');
    set(handles.runbtn2,'Enable','on');
    set(handles.runbothbtn,'Enable','on');
end
