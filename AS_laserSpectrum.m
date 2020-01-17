% zeke barger 011720 
% plots EEG power spectrum to laser onset
% simply averages across all trials
% takes a dirlist or mouselist... but a dirlist is probably better
% this re-uses some code from AS_brainstate so it's a little more
% complicated than it needs to be :/

function [] = AS_laserSpectrum()
%% set these parameters before you run this
SR = 512; % EEG sampling rate (e.g., 512)
before = 240; % seconds before laser onset to display (default 240)
after = 240; % seconds after to display (default 240)
laser_duration = 120; % laser duration in sec
spectrogram_window = 5; % spectrogram window, in seconds (e.g., 5)
spectrogram_step = 2.5; % spectrogram time step (e.g., 2.5)
spectrogram_freq_lo = 1; % min frequency to show
spectrogram_freq_hi = 30; % max frequency to show
z_min = -4; % lower limit on z-scored values in the plot
z_max = 4; % upper limit

%% load and process data
% EEG samples to collect before, during and after laser
before_samples = round(before*SR); % number before laser
after_samples = round(after*SR); % number after laser
laser_samples = round(laser_duration*SR); % number during laser

% have user select the mouselist 
[fname,fpath,~] = uigetfile('*','Select mouselist or dirlist');
% check if something was selected
if ~ischar(fname)
    disp('no list selected');
    return
end

% Load either a mouselist (for multiple mice) or dirlist (for 1 mouse).
% If only a dirlist is selected, the analysis will be the same except we
% will will treat each day's data as if it came from a different mouse.
% This lets us re-use the same code. dirlistMode indicates whether we are
% dealing with a dirlist (1) or a mouselist (0).
temp = load([fpath,fname],'mouselist'); % load list of dirlists
if ~isfield(temp,'mouselist')
    temp = load([fpath,fname],'dirlist'); % load dirlist
    if ~isfield(temp,'dirlist')
        return
    end
    mouselist = temp.dirlist;
    dirlistMode = 1;
else
    mouselist = temp.mouselist;
    dirlistMode = 0;
end

if isempty(mouselist) % nothing in the list
    disp('the list is empty');
    return
end

% load data
nMice = length(mouselist); % how many mice there are
allMiceData = cell(1, nMice); % laser locked brain states, one cell per mouse
allMiceFiles = cell(1, nMice); % locations of files
nDays = zeros(1,nMice); % how many days there are for each mouse
% if we are dealing with a mouselist
if ~dirlistMode
    % check that the contents of the mouselist are ok, and gather filenames
    for m = 1:nMice % for each mouse
        % find EEG and laser data
        temp = load(mouselist{m},'dirlist');
        dirlist = temp.dirlist;
        [problemString, fileNames] = AS_checkList(dirlist, {'EEG','laser'});
        % if something is wrong with this mouse
        if ~isempty(problemString)
            % display the problem and quit
            disp(problemString)
            return
        end
        % continue on, storing the locations of the files for this mouse
        allMiceFiles{m} = fileNames;
        nDays(m) = length(dirlist);
    end
else % we are dealing with a dirlist
    [problemString, fileNames] = AS_checkList(mouselist, {'EEG','laser'});
    % if something is wrong with this list
    if ~isempty(problemString)
        % display the problem and quit
        disp(problemString)
        return
    end
    for m = 1:nMice % for each recording
        allMiceFiles{m} = {fileNames{m}};
        nDays(m) = 1; % just one "day"
    end
end

trialsPerMouse = zeros(1,nMice);
wb = waitbar(0,'Loading...');
for m = 1:nMice % for each mouse
    waitbar(m/nMice,wb); % update progress bar
    allMiceData{m} = []; % initialize array of trial info
    for d = 1:nDays(m) % for each day
        % get laser-locked EEG
        allMiceData{m} = [allMiceData{m}; getLaserEEG(allMiceFiles{m}{d}{1}, ...
            allMiceFiles{m}{d}{2},before_samples,laser_samples,after_samples, SR)];
    end
    % number of rows is the number of trials for this animal
    trialsPerMouse(m) = size(allMiceData{m}, 1);
end
delete(wb)

% no trials found? quit
if sum(trialsPerMouse) == 0
    disp('no trials found :(')
    return
end

% preallocate a huge matrix of spectrograms
% by first finding the size of the spectrogram
% but first, find first trial and get its EEG
EEG = allMiceData{find(trialsPerMouse, 1)}(1,:);
[s_temp, ~, ~] = makeSpectro(EEG, SR, spectrogram_step, ...
    spectrogram_window, spectrogram_freq_lo, spectrogram_freq_hi);
specs = zeros(size(s_temp,1), size(s_temp,2), sum(trialsPerMouse));

% for each mouse, get spectrograms
idx = 1; % nth spectrogram created
wb = waitbar(0,'creating spectrograms...');
for i = 1:nMice % mouse
    for j = 1:trialsPerMouse(i) % trial
        waitbar(idx/sum(trialsPerMouse),wb); % update progress bar
        [specs(:,:,idx), t, f] = makeSpectro(allMiceData{i}(j,:), SR, ...
            spectrogram_step, spectrogram_window, spectrogram_freq_lo,...
            spectrogram_freq_hi);
        idx = idx+1;
    end
end
delete(wb)

% determine when the laser starts in the spectrograms
percent_before = before_samples / (before_samples + after_samples + laser_samples);
t_steps_before = round(percent_before * length(t));

% for each frequency band, z-score based on pre-laser average / SD
mu = zeros(length(f),1);
sig = zeros(1,length(f));
for i = 1:length(f)
    mu(i) = mean(mean(specs(i,1:t_steps_before,:))); 
    sig(i) = std(reshape(specs(i,1:t_steps_before,:), 1, ...
        numel(specs(i,1:t_steps_before,:))));
end
spec_avg = mean(specs,3);
spec_z = spec_avg - mu;
spec_z = spec_z ./ sig;
spec_z(spec_z < z_min) = z_min;
spec_z(spec_z > z_max) = z_max;

%% plot
figure('Color','w');
imagesc(t,f,spec_z);




% Make a matrix (M) of EEG data (rows are trials, columns are
%   timepoints) time locked to laser stimuli
% Trials that are too short are discarded
% inputs:
% EEGfileName: points to the EEG file
% laserFileName: points to the laser file
% <something>Samples: number of samples before, during, and after each laser pulse
% SR: EEG sampling rate
function M = getLaserEEG(EEGFileName, laserFileName, beforeSamples, ...
    laserSamples, afterSamples, SR)
% load EEG and laser onsets/offsets
temp = load(EEGFileName,'EEG');
EEG = temp.EEG;
temp = load(laserFileName,'laser');
laser = temp.laser;

% find sample # of laser onsets
idx_start = floor(laser(:,1) * SR) + 1;
% number of stimuli
ntrials = length(idx_start); 

% matrix of EEG data on each trial
M = nan(ntrials,laserSamples+beforeSamples+afterSamples); 

% collect data
i = 1;
while i <= size(M,1)
    % if certain conditions aren't met, discard trial
    if (~(idx_start(i)-beforeSamples > 0) ||... % trial begins before start of recording
            ~(idx_start(i)+laserSamples+afterSamples-1 <...
            length(EEG))) % or ends after the end
        M(i,:) = []; % delete row in M
        idx_start(i) = []; % delete laser time
        ntrials = ntrials-1; % update trial count
        continue % don't update i
    else % store the data from the trial
        M(i,:) = EEG(idx_start(i)-beforeSamples:idx_start(i)+...
            laserSamples+afterSamples-1);
        i = i+1;
    end
end

% create a multi-taper spectrogram of EEG data
% s is the spectrogram, t is the time axis, f is the frequency axis
% EEG is the EEG data, SR its sampling rate, winstep is window step,
% winwid is window width
function [s, t, f] = makeSpectro(EEG, SR, winstep, winwid, lo, hi)
params = struct;
params.pad = -1;
params.Fs = SR;
params.fpass = [lo hi];
params.tapers = [3 5];
% make sure EEG is a row
if ~isrow(EEG)
    EEG = EEG';
end
% truncate EEG to a multiple of SR*winstep
EEG = EEG(1:(length(EEG)-mod(length(EEG), SR*winstep)));
% pad the EEG signal so that the first bin starts at time 0
EEG = [EEG(1:round(SR*(winwid-winstep)/2)), EEG, ...
    EEG((end+1-round(SR*(winwid-winstep)/2)):end)];
[s, t, f] = mtspecgramc(EEG, [winwid, winstep], params);
% adjust time axis to reflect this change
t = t - (winwid-winstep)/2;