% zeke barger 081120
% Convert data from Min's format into something AccuSleep can use
% inputs:
% min_file: name of a file in Min's format
% output_folder: name of the folder where converted files will be saved
% note: you should use a new folder for every Min-style file
% also note: this doesn't get EEG and EMG data, so it's not possible to
% manually re-annotate these outputs

function AS_convertMin(min_file, output_folder)
% load the data
M = load(min_file);

% extract brain state labels
labels = M.bCluster.rem*1 + M.bCluster.wake*2 + M.bCluster.sws*3;

% extract laser timing
laser_start_idx = [1; find(diff(M.recordingFile.LasT(:,1))>10)+1];
laser_end_idx = [find(diff(M.recordingFile.LasT(:,1))>10); size(M.recordingFile.LasT,1)];
laser = [M.recordingFile.LasT(laser_start_idx,1),...
    M.recordingFile.LasT(laser_end_idx,1)];

% save these
save([output_folder,filesep,'labels.mat'],'labels');
save([output_folder,filesep,'laser.mat'],'laser');
