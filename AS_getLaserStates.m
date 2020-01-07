% zeke barger 121319
% Make a matrix (M) of brain state data (rows are trials, columns are
%   timepoints, entries are brain states) time locked to laser stimuli
% Trials that are too short or contain undefined epochs are discarded
% inputs:
% fileName: points to the labels file
% laserFileName: points to the laser file
% <something>bins: number of epochs before, during, and after each laser pulse
% dt: length of epochs for brain states

function M = AS_getLaserStates(labelFileName, laserFileName, beforeBins, laserBins, afterBins, dt)
% load brain state labels and laser onsets/offsets
temp = load(labelFileName,'labels');
labels = temp.labels;
temp = load(laserFileName,'laser');
laser = temp.laser;

% convert any undefined epochs to NaN
labels(labels > 3 | labels < 1) = NaN;

% find epoch # of laser onsets
idx_start = floor(laser(:,1) / dt) + 1;
% number of stimuli
ntrials = length(idx_start); 

% matrix of brain states on each trial
M = nan(ntrials,laserBins+beforeBins+afterBins); 

% collect data
i = 1;
while i <= size(M,1)
    % if certain conditions aren't met, discard trial
    if (~(idx_start(i)-beforeBins > 0) ||... % trial begins before start of recording
            ~(idx_start(i)+laserBins+afterBins-1 < length(labels))) ||... % or ends after the end
            (any(isnan((labels(idx_start(i)-beforeBins:idx_start(i)+... % or contains NaNs
            laserBins+afterBins-1)))))
        M(i,:) = []; % delete row in M
        idx_start(i) = []; % delete laser time
        ntrials = ntrials-1; % update trial count
        continue % don't update i
    else % store the data from the trial
        M(i,:) = labels(idx_start(i)-beforeBins:idx_start(i)+...
            laserBins+afterBins-1);
        i = i+1;
    end
end