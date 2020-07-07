function varargout = AS_exportIntan(varargin)
% zeke barger 070720
% Exports EEG/EMG recording from Intan system into a format that can be
% interpreted by AccuSleep. The onsets and offsets of laser stimuli, if
% present are saved (in seconds).

% TODO: don't overwrite folders when making them
% partial files

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AS_exportIntan_OpeningFcn, ...
    'gui_OutputFcn',  @AS_exportIntan_OutputFcn, ...
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


% --- Executes just before AS_exportIntan is made visible.
function AS_exportIntan_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for AS_exportIntan
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
setappdata(handles.C,'laser_channel','ADC1');

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
set(handles.adcBox,'String','ADC1');

function varargout = AS_exportIntan_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function blockbtn_Callback(~, ~, handles)
% choose the block
blockPath = uigetdir(get(handles.blockbox,'String'), 'Choose Intan folder to process');

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
drawnow;

for n = a:b % for each mouse selected
    disptext(handles, ['Processing mouse ',num2str(n),', please wait...']);
    drawnow;
    
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
    if ~strcmp(location(end),'\')
        location = [location,'\'];
    end
    % opposite for output folder name
    if contains(outputName,'\')
        errordlg(['Output folder name should not contain slashes. ',...
            'A good folder name might be the date of the experiment']);
        setListLock(handles,0);
        return
    end
    
    % load the data
    % First, get settings file
    all_files = struct2cell(dir(blockPath))';
    isfolder = cell2mat(all_files(:,5)); % find entries that are folders
    all_files = all_files(~isfolder,1); % remove those
    if any(contains(all_files,'.xml'))
        settingsFilename = all_files{contains(all_files,'xml')};
    else
        errordlg('Could not find a .xml file containing recording settings.');
        setListLock(handles,0);
        return
    end
    
    originalSR = getSR([blockPath,'\',settingsFilename]); % get sampling rate
    
    % try to load the data
    eegMatch = contains(all_files,[string(['CH',num2str(EEGch),'.']), ...
        string(['CH',num2str(EEGch),'_'])]);
    if any(eegMatch)
        eegFilename = all_files{eegMatch};
    else
        errordlg('Could not find the EEG file.');
        setListLock(handles,0);
        return
    end
    emgMatch = contains(all_files,[string(['CH',num2str(EMGch),'.']), ...
        string(['CH',num2str(EMGch),'_'])]);
    if any(emgMatch)
        emgFilename = all_files{emgMatch};
    else
        errordlg('Could not find the EMG file.');
        setListLock(handles,0);
        return
    end
    laser_channel = getappdata(handles.C,'laser_channel');
    laserMatch = contains(all_files,[string([laser_channel,'.']), ...
        string([laser_channel,'_'])]);
    if any(laserMatch)
        laserFilename = all_files{laserMatch};
    else
        errordlg('Could not find the laser file.');
        setListLock(handles,0);
        return
    end
    
    % actually load them all
    try
        [eegData, ~, ~] = load_open_ephys_data([blockPath,'\',eegFilename]);
        [emgData, ~, ~] = load_open_ephys_data([blockPath,'\',emgFilename]);
        [laserData, ~, ~] = load_open_ephys_data([blockPath,'\',laserFilename]);
    catch
        errordlg('Error loading files.');
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
    
    % downsample EEG, EMG
    EEG = eegData(1:(originalSR/SR):end);
    EMG = emgData(1:(originalSR/SR):end);
    
    % get laser stimulus onsets and offsets
    [onsets, offsets] = binarizeLaser(laserData, originalSR);
    
    % new way of saving laser data (more compact)
    if ~isempty(onsets) && ~isempty(offsets)
        laser = [onsets, offsets];
        save([location,outputName,'\laser.mat'],'laser');
    end
    
    save([location,outputName,'\EEG.mat'],'EEG');
    save([location,outputName,'\EMG.mat'],'EMG');
    
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
    set(handles.prevBtn,'Enable','off');
else
    set(handles.blockbtn,'Enable','on');
    set(handles.blockbox,'Enable','on');
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
    set(handles.prevBtn,'Enable','on');
end

% get sampling rate from a settings xml file
function sr = getSR(fname)
a = xmlread(fname);
c = a.getElementsByTagName('EDITOR');
d = c.item(0);
e = d.getAttributes;
srString=[];
for i = 0:d.getAttributes.getLength
    if ~isempty(strfind(e.item(i), 'SampleRate'))
        srString = char(e.item(i));
        break
    end
end
qidx = find(srString == '"');
sr = srString(qidx(1)+1 : qidx(2)-1);
sr = str2num(sr); % this is actually the index in the following samplerate list
possibilities = [1 1.25 1.5 2 2.5 3 3.33 4 5 6.25 8 10 12.5 15 20 25 30];
sr = round(possibilities(sr) * 1000);

function [onsets, offsets] = binarizeLaser(laser_data, SR)
thresh = .98; % noise threshold compared to max value
minISI = 30; % max allowable interval between pulses in a train, in seconds
minSamplesInPulse = 4;

m = max(laser_data);
laser_data(laser_data < thresh*m) = 0; % remove low values
laser_data(laser_data > 0) = 1; % set high values to 1
f = find(laser_data);
if isempty(f)
    onsets = [];
    offsets = [];
    return
end
% remove singletons
laser_data(f(1)) = 0;
laser_data(f(end)) = 0;
for i = 2:(length(f)-1)
    if laser_data(f(i)-1) == 0 && laser_data(f(i)+1) == 0
        laser_data(f(i)) = 0;
    end
end

minISI_in_samples = minISI*SR;

d = diff(laser_data);
pulse_onsets = find(d==1);
pulse_offsets = find(d==-1);
% check if first offset is before first onset, and remove if so
while pulse_offsets(1) < pulse_onsets(1)
    pulse_offsets(1) = [];
end
% check if last onset is after last offset, and remove if so
while pulse_onsets(end) > pulse_offsets(end)
    pulse_onsets(end) = [];
end
onsets = [pulse_onsets(1);pulse_onsets(1+find(diff(pulse_onsets) > minISI_in_samples))];
offsets = [pulse_offsets(find(diff(pulse_offsets) > minISI_in_samples)); pulse_offsets(end)];

% remove pairs that are too close together
i = 1;
while i <= length(onsets)
    if offsets(i) - onsets(i) < minSamplesInPulse
        offsets(i) = [];
        onsets(i) = [];
    else
        i = i + 1;
    end
end

onsets = (onsets+1)/SR;
offsets = (offsets+1)/SR;

function adcBox_Callback(hObject, eventdata, handles)
% update the stored laser channel
setappdata(handles.C,'laser_channel',get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function adcBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in prevBtn.
function prevBtn_Callback(hObject, eventdata, handles)
% determine if path is valid
blockPath = getappdata(handles.C,'blockPath');
if ~exist(blockPath,'dir')  
    ('No folder selected');
    return
end

% load the data, if it's not already loaded
setListLock(handles,1);

% First, get settings file
all_files = struct2cell(dir(blockPath))';
isfolder = cell2mat(all_files(:,5)); % find entries that are folders
all_files = all_files(~isfolder,1); % remove those
if any(contains(all_files,'.xml'))
    settingsFilename = all_files{contains(all_files,'xml')};
else
    errordlg('Could not find a .xml file containing recording settings.');
    setListLock(handles,0);
    return
end

originalSR = getSR([blockPath,'\',settingsFilename]); % get sampling rate

% try to load the data
chMatch = contains(all_files,"CH");
if any(chMatch)
    chFilenames = all_files(chMatch);
else
    errordlg('Could not find any data.');
    setListLock(handles,0);
    return
end

% actually load them all
chdata = {};
try
    for i = 1:length(chFilenames)
        chdata{i} = load_open_ephys_data([blockPath,'\',chFilenames{i}],...
            'Indices',1:(originalSR*60*20));
    end
catch
    errordlg('Error loading files.');
    setListLock(handles,0);
    return
end
drawnow;

durationInSeconds = length(chdata{1})/originalSR;

disptext(handles, 'Finished loading.')
drawnow;

% get information about the data
numChannels = length(chdata);
numSamples = length(chdata{1});
% calculate min, max, text position for each channel
stats = zeros(numChannels, 3);
for i = 1:numChannels
    stats(i,1) = min(chdata{i});
    stats(i,2) = max(chdata{i});
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
plot(chdata{1}(1:(originalSR/200):end));
labels{1} = chFilenames{1};
for i = 2:numChannels
    hold on, plot(chdata{i}(1:(originalSR/200):end)-...
        (stats(i,2)-stats(i-1,1)));
    stats(i,:) = stats(i,:) - (stats(i,2)-stats(i-1,1));
    labels{i} = chFilenames{i};
end
axis tight
title('Preview (~20 mins)')
set(gca,'YTick',flipud(stats(:,3)),'YTickLabel',fliplr(labels),'XTick',[])
set(gca,'TickLabelInterpreter','none')
warning('on','MATLAB:colon:nonIntegerIndex');
setListLock(handles,0);


function [data, timestamps, info] = load_open_ephys_data(filename,varargin)
%     Copyright (C) 2014 Open Ephys
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.


p = inputParser;
p.addRequired('filename');
p.addParameter('Indices',[],@(x) isempty(x) || isrow(x) && all(diff(x) > 0) ...
    && all(fix(x) == x) && all(x > 0));
p.parse(filename,varargin{:});

range_pts = p.Results.Indices;

filetype = filename(max(strfind(filename,'.'))+1:end); % parse filetype

fid = fopen(filename);
filesize = getfilesize(fid);

% constants
NUM_HEADER_BYTES = 1024;
SAMPLES_PER_RECORD = 1024;
RECORD_MARKER = [0 1 2 3 4 5 6 7 8 255]';
RECORD_MARKER_V0 = [0 0 0 0 0 0 0 0 0 255]';

% constants for pre-allocating matrices:
MAX_NUMBER_OF_SPIKES = 1e6;
MAX_NUMBER_OF_RECORDS = 1e6;
MAX_NUMBER_OF_CONTINUOUS_SAMPLES = 1e8;
MAX_NUMBER_OF_EVENTS = 1e6;
SPIKE_PREALLOC_INTERVAL = 1e6;

if strcmp(filetype, 'continuous')
    % https://open-ephys.atlassian.net/wiki/spaces/OEW/pages/65667092/Open+Ephys+format#OpenEphysformat-Continuousdatafiles(.continuous)
    disp(['Loading ' filename '...']);
    index = 0;
    hdr = fread(fid, NUM_HEADER_BYTES, 'char*1');
    eval(char(hdr'));
    info.header = header;
    if (isfield(info.header, 'version'))
        version = info.header.version;
    else
        version = 0.0;
    end
    % pre-allocate space for continuous data
    data = zeros(MAX_NUMBER_OF_CONTINUOUS_SAMPLES, 1);
    info.ts = zeros(1, MAX_NUMBER_OF_RECORDS);
    info.nsamples = zeros(1, MAX_NUMBER_OF_RECORDS);
    
    if version >= 0.2
        info.recNum = zeros(1, MAX_NUMBER_OF_RECORDS);
    end
    current_sample = 0;
    RECORD_SIZE = 10 + SAMPLES_PER_RECORD*2 + 10; % size of each continuous record in bytes
    if version >= 0.2
        RECORD_SIZE = RECORD_SIZE + 2; % include recNum
    end
    if isempty(range_pts)
        while ftell(fid) + RECORD_SIZE <= filesize % at least one record remains
            go_back_to_start_of_loop = 0;
            index = index + 1;
            if (version >= 0.1)
                timestamp = fread(fid, 1, 'int64', 0, 'l');
                nsamples = fread(fid, 1, 'uint16',0,'l');
                if version >= 0.2
                    recNum = fread(fid, 1, 'uint16');
                end
            else
                timestamp = fread(fid, 1, 'uint64', 0, 'l');
                nsamples = fread(fid, 1, 'int16',0,'l');
            end
            if nsamples ~= SAMPLES_PER_RECORD && version >= 0.1
                disp(['  Found corrupted record...searching for record marker.']);
                % switch to searching for record markers
                last_ten_bytes = zeros(size(RECORD_MARKER));
                for bytenum = 1:RECORD_SIZE*5
                    byte = fread(fid, 1, 'uint8');
                    last_ten_bytes = circshift(last_ten_bytes,-1);
                    last_ten_bytes(10) = double(byte);
                    if last_ten_bytes(10) == RECORD_MARKER(end)
                        sq_err = sum((last_ten_bytes - RECORD_MARKER).^2);
                        if (sq_err == 0)
                            disp(['   Found a record marker after ' int2str(bytenum) ' bytes!']);
                            go_back_to_start_of_loop = 1;
                            break; % from 'for' loop
                        end
                    end
                end
                
                % if we made it through the approximate length of 5 records without
                % finding a marker, abandon ship.
                if bytenum == RECORD_SIZE*5
                    disp(['Loading failed at block number ' int2str(index) '. Found ' ...
                        int2str(nsamples) ' samples.'])
                    break; % from 'while' loop
                end 
            end
            if ~go_back_to_start_of_loop
                block = fread(fid, nsamples, 'int16', 0, 'b'); % read in data
                fread(fid, 10, 'char*1'); % read in record marker and discard
                data(current_sample+1:current_sample+nsamples) = block;
                current_sample = current_sample + nsamples;
                info.ts(index) = timestamp;
                info.nsamples(index) = nsamples;
                if version >= 0.2
                    info.recNum(index) = recNum;
                end
            end
        end
    elseif ~isempty(range_pts)
        if (version >= 0.1)
            if version >= 0.2
                m = memmapfile(filename,....
                    'Format',{'int64',1,'timestamp';...
                    'uint16',1,'nsamples';...
                    'uint16',1,'recNum';...
                    'int16',[SAMPLES_PER_RECORD, 1],'block';...
                    'uint8',[1, 10],'marker'},...
                    'Offset',NUM_HEADER_BYTES,'Repeat',Inf);
            else
                m = memmapfile(filename,....
                    'Format',{'int64',1,'timestamp';...
                    'uint16',1,'nsamples';...
                    'int16',[SAMPLES_PER_RECORD, 1],'block';...
                    'uint8',[1, 10],'marker'},...
                    'Offset',NUM_HEADER_BYTES,'Repeat',Inf); %TODO not tested
            end
        else
            m = memmapfile(filename,....
                'Format',{'uint64',1,'timestamp';...
                'int16',1,'nsamples';...
                'int16',[SAMPLES_PER_RECORD, 1],'block';...
                'uint8',[1, 10],'marker'},...
                'Offset',NUM_HEADER_BYTES,'Repeat',Inf); %TODO not tested
        end
        tf = false(length(m.Data)*SAMPLES_PER_RECORD,1);
        tf(range_pts) = true;
        Cblk = cell(length(m.Data),1);
        Cts = cell(length(m.Data),1);
        Cns = cell(length(m.Data),1);
        Ctsinterp = cell(length(m.Data),1);
        if version >= 0.2
            Crn = cell(length(m.Data),1);
        end
        for i = 1:length(m.Data)
            if any(tf(SAMPLES_PER_RECORD*(i-1)+1:SAMPLES_PER_RECORD*i))
                Cblk{i} = m.Data(i).block(tf(SAMPLES_PER_RECORD*(i-1)+1:SAMPLES_PER_RECORD*i));
                Cts{i} = double(m.Data(i).timestamp);
                Cns{i} = double(m.Data(i).nsamples);
                if version >= 0.2
                    Crn{i} = double(m.Data(i).recNum);
                end
                tsvec = Cts{i}:Cts{i}+Cns{i}-1;
                Ctsinterp{i} = tsvec(tf(SAMPLES_PER_RECORD*(i-1)+1:SAMPLES_PER_RECORD*i));
            end
        end
        data = vertcat(Cblk{:});
        data = double(swapbytes(data)); % big endian
        info.ts = [Cts{:}];
        info.nsamples = [Cns{:}];
        if version >= 0.2
            info.recNum = [Crn{:}];
        end
        index = length(info.nsamples);
        current_sample = length(range_pts);
        timestamps = [Ctsinterp{:}]';
        %TODO check for corrupted file, not tested
        if any(info.nsamples ~= SAMPLES_PER_RECORD) && version >= 0.1
            disp(['  Found corrupted record...searching for record marker.']);
            k = find(info.nsamples ~= SAMPLES_PER_RECORD,1,'first');
            ns = info.nsamples(k);
            if version >= 0.2
                offset = NUM_HEADER_BYTES + RECORD_SIZE * (k-1) + 12;
            else
                offset = NUM_HEADER_BYTES + RECORD_SIZE * (k-1) + 10;
            end
            status = fseek(fid,offset,'bof');
            last_ten_bytes = zeros(size(RECORD_MARKER));
            for bytenum = 1:RECORD_SIZE*5
                byte = fread(fid, 1, 'uint8');
                last_ten_bytes = circshift(last_ten_bytes,-1);
                last_ten_bytes(10) = double(byte);
                if last_ten_bytes(10) == RECORD_MARKER(end)
                    sq_err = sum((last_ten_bytes - RECORD_MARKER).^2);
                    if (sq_err == 0)
                        disp(['   Found a record marker after ' int2str(bytenum) ' bytes!']);
                        go_back_to_start_of_loop = 1;
                        break; % from 'for' loop
                    end
                end
            end
            if bytenum == RECORD_SIZE*5
                disp(['Loading failed at block number ' int2str(index) '. Found ' ...
                    int2str(nsamples) ' samples.'])
            end
        end
    end
    % crop data to the correct size
    data = data(1:current_sample);
    info.ts = info.ts(1:index);
    info.nsamples = info.nsamples(1:index);
    if version >= 0.2
        info.recNum = info.recNum(1:index);
    end
    % convert to microvolts
    data = data.*info.header.bitVolts;
    if isempty(range_pts)
        timestamps = nan(size(data));
        current_sample = 0;
        if version >= 0.1
            for record = 1:length(info.ts)
                ts_interp = info.ts(record):info.ts(record)+info.nsamples(record);
                timestamps(current_sample+1:current_sample+info.nsamples(record)) = ts_interp(1:end-1);
                current_sample = current_sample + info.nsamples(record);
            end
        else % v0.0; NOTE: the timestamps for the last record will not be interpolated
            for record = 1:length(info.ts)-1
                ts_interp = linspace(info.ts(record), info.ts(record+1), info.nsamples(record)+1);
                timestamps(current_sample+1:current_sample+info.nsamples(record)) = ts_interp(1:end-1);
                current_sample = current_sample + info.nsamples(record);
            end
        end
    end
else
    error('File extension not recognized. Please use a ''.continuous'', ''.spikes'', or ''.events'' file.');    
end

fclose(fid); % close the file
if (isfield(info.header,'sampleRate'))
    if ~ischar(info.header.sampleRate)
        timestamps = timestamps./info.header.sampleRate; % convert to seconds
    end
end

function filesize = getfilesize(fid)
fseek(fid,0,'eof');
filesize = ftell(fid);
fseek(fid,0,'bof');