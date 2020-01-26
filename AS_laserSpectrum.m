% zeke barger 012620 
% plots EEG power spectrum locked to laser stimuli
% takes a dirlist or mouselist
% for each recording (day), find z-scored changes from baseline
% for each mouse, average across days weighted by trials per day
% then average across mice

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
z_min = -1; % lower limit on z-scored values in the plot (e.g., -1)
z_max = 1; % upper limit

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
temp = load([fpath,fname],'mouselist'); % load list of dirlists
if ~isfield(temp,'mouselist')
    temp = load([fpath,fname],'dirlist'); % load dirlist
    if ~isfield(temp,'dirlist')
        return
    end
    mouselist = {[fpath,fname]};
else
    mouselist = temp.mouselist;
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

trialsPerDay = zeros(nMice, max(nDays));
wb = waitbar(0,'Loading...');
for m = 1:nMice % for each mouse
    waitbar(m/nMice,wb); % update progress bar
    allMiceData{m} = {}; % initialize 
    for d = 1:nDays(m) % for each day
        % get laser-locked EEG
        allMiceData{m}{d} = getLaserEEG(allMiceFiles{m}{d}{1}, ...
            allMiceFiles{m}{d}{2},before_samples,laser_samples,after_samples, SR);
        % number of rows is the number of trials for this animal
        trialsPerDay(m,d) = size(allMiceData{m}{d}, 1);
    end
end
delete(wb)

% no trials found? quit
totalNumTrials = sum(sum(trialsPerDay));
if totalNumTrials == 0
    disp('no trials found :(')
    return
end

% preallocate spectrograms by first finding the size of the spectrogram
% so, find first trial and get its EEG
[r,c] = find(trialsPerDay); % select days with trials
EEG = allMiceData{r(1)}{c(1)}(1,:); % get EEG of first trial
[s_temp, ~, ~] = makeSpectro(EEG, SR, spectrogram_step, ...
    spectrogram_window, spectrogram_freq_lo, spectrogram_freq_hi);
specs = cell(1,nMice);
wb = waitbar(0,'creating spectrograms...');
idx = 1; % nth spectrogram created

for m = 1:nMice
    for d = 1:nDays(m) % for each day, preallocate
        specs{m}{d} = zeros(size(s_temp,1), size(s_temp,2), trialsPerDay(m,d));
        for j = 1:trialsPerDay(m,d) % trial
            waitbar(idx/totalNumTrials,wb); % update progress bar
            [specs{m}{d}(:,:,j), t, f] = makeSpectro(allMiceData{m}{d}(j,:), SR, ...
                spectrogram_step, spectrogram_window, spectrogram_freq_lo,...
                spectrogram_freq_hi);
            idx = idx+1;
        end
    end
end
delete(wb)

% determine when the laser starts in the spectrograms
percent_before = before_samples / (before_samples + after_samples + laser_samples);
t_steps_before = round(percent_before * length(t));

% preallocate z-scored spectrogram stacks for each mouse
zStacks = cell(1, nMice);
mouseMeans = zeros(size(s_temp,1), size(s_temp,2), nMice);
% for each frequency band, z-score based on pre-laser average / SD
for m = 1:nMice
    zStacks{m} = zeros(size(s_temp,1), size(s_temp,2), nDays(m));
    for d = 1:nDays(m)
        mu = zeros(1,length(f));
        sig = zeros(1, length(f));
        for i = 1:length(f)
            mu(i) = mean(mean(specs{m}{d}(1:t_steps_before,i,:)));
            sig(i) = std(reshape(specs{m}{d}(1:t_steps_before,i,:), 1, ...
                numel(specs{m}{d}(1:t_steps_before,i,:))));
        end
        spec_avg = mean(specs{m}{d},3);
        spec_z = spec_avg - mu;
        spec_z = spec_z ./ sig;
        % for weighted averaging
        spec_z = spec_z * trialsPerDay(m,d);
        zStacks{m}(:,:,d) = spec_z;
    end
    mouseMeans(:,:,m) = mean(sum(zStacks{m},3)/sum(trialsPerDay(m,:)),3);
end
grandMean = mean(mouseMeans,3);
%% plot
t = t - before;
figure('Color','w');
imagesc(t,f,grandMean');
axis xy
caxis([z_min z_max])
colormap(gca,cat(1,cat(2,(0:127)'./127,(0:127)'./127,ones(128,1)),...
    cat(2,ones(128,1),(127:-1:0)'./127,(127:-1:0)'./127)))
ylabel('Frequency (Hz)')
xlabel('Time (s)')
colorbar;
title('z-scored EEG spectrogram')
yl = ylim;
hold on, plot([0 0],yl,'k');
hold on, plot([laser_duration,laser_duration],yl,'k');

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
