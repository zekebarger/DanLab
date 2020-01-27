% zeke barger 012720
% plots brain state relative to laser onset
% averages across animals and uses bootstrapping
% requires a mouselist (list of dirlist file locations)

%TODO
% transition asterisks and diagram
% eyfp comparison


function [] = AS_brainstate()
%% set parameters
iters = 1000; % iterations for generating bootstrap confidence intervals (default 10000)
before = 240; % seconds before laser onset to display (default 240)
after = 240; % seconds after to display (default 240)
laser_duration = 120; % laser duration in sec
state_bin = 2.5; % epoch length for brain state classification
% order to plot the lines (last entry is plotted on top). 1=R, 2=W, 3=N
plot_order = [2 3 1];
% line colors (REM, wake, NREM)
colors = [43 160 220;
    129 131 132
    255 194 61] ./ 255;
% error bar colors
errorcolors = [43 160 220;
    81 82 82
    255 194 61] ./ 255;
% color for edge of error bars
edgecolors = [150 207 238;
    194 195 196;
    244 212 153] ./ 255;
% transparency of error bars
bar_alphas = [.4 .4 .4];


% no need to edit these lines
% epochs to plot before, during and after laser
before_epochs = round(before/state_bin); % number of bins before laser
after_epochs = round(after/state_bin); % number of bins after laser
laser_epochs = round(laser_duration / state_bin); % number of bins during laser


%% load and process data
% have user select the mouselist
[fname,fpath,~] = uigetfile('*','Select mouselist or dirlist for TREATMENT condition');
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
        % find sleep stage and laser data
        temp = load(mouselist{m},'dirlist');
        dirlist = temp.dirlist;
        [problemString, fileNames] = AS_checkList(dirlist, {'labels','laser'});
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
    [problemString, fileNames] = AS_checkList(mouselist, {'labels','laser'});
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

wb = waitbar(0,'Loading...');
for m = 1:nMice % for each mouse
    waitbar(m/nMice,wb); % update progress bar
    allMiceData{m} = []; % initialize array of trial info
    for d = 1:nDays(m) % for each day
        % get laser-locked brain states
        allMiceData{m} = [allMiceData{m}; AS_getLaserStates(allMiceFiles{m}{d}{1}, ...
            allMiceFiles{m}{d}{2},before_epochs,laser_epochs,after_epochs, state_bin)];
    end
end
delete(wb)

% for each mouse, average across days
mouseAvgs = []; % state, time, mouse
for i = 1:3 % state
    for j = 1:nMice % mouse
        mouseAvgs(i,:,j) = mean(allMiceData{j}==i);
    end
end

%% find mean and CIs
trialsPerMouse = zeros(1,nMice); % number of trials for each mouse
for i = 1:nMice
    trialsPerMouse(i) = size(allMiceData{i},1);
end
% totaltrials = sum(mousetrials); % total number of trials in the whole experiment
trialdur = size(allMiceData{1},2); % number of time bins in each trial
bsavg = zeros(iters,trialdur,3); % iter by time by state bootstrapped likelihood

wb = waitbar(0,'Bootstrapping...');
diffs = zeros(iters, 3); % pre vs laser difference
pre_lsr_pts = before / state_bin;
lsr_pts = laser_duration / state_bin;
for itr = 1:iters
    waitbar(itr/iters,wb);
    bsItrAvg = zeros(nMice,trialdur,3);
    for m = 1:nMice % get data randomly from each mouse
        bsdata = allMiceData{m}(randi(trialsPerMouse(m),1,trialsPerMouse(m)),:);
        % average across trials for this mouse
        for k = 1:3 % state
            bsItrAvg(m,:,k) = mean(bsdata==k);
        end
    end
    
    bsavg(itr,:,:) = mean(bsItrAvg,1);
    diffs(itr,:) = squeeze(mean(bsavg(itr,1:pre_lsr_pts,:))) -...
        squeeze(mean(bsavg(itr,(1:lsr_pts)+pre_lsr_pts,:)));
    
end
delete(wb);

CIs = prctile(bsavg,[2.5 97.5],1);
CIs = CIs*100;

grandMean = 100*mean(mouseAvgs,3);

%% plot
state_labels = {'REM','Wake','NREM'};
state_labels = state_labels(plot_order);

figure('Color','w')
fill([0 laser_duration laser_duration 0],[0 0 100 100],[154 154 255]./255,'EdgeAlpha',0,'FaceAlpha',.5);

t = (1:(after_epochs+before_epochs+laser_epochs))*state_bin - before_epochs*state_bin - state_bin/2;
t_flip = cat(2,t,fliplr(t));

legend_mean = [];
legend_shading = [];
legend_top = [];
legend_bottom = [];
for i = plot_order
    hold on, legend_shading(i) = fill(t_flip, cat(2,CIs(1,:,i),fliplr(CIs(2,:,i))),...
        errorcolors(i,:), 'EdgeAlpha', 0,'FaceAlpha', bar_alphas(i),'LineWidth',1.5);
    hold on, legend_mean(i) = plot(t, grandMean(i,:), 'Color',colors(i,:), 'LineWidth',1.5);
    hold on, legend_top(i) = plot(t, CIs(2,:,i), 'Color',edgecolors(i,:), 'LineWidth',1);
    hold on, legend_bottom(i) = plot(t, CIs(1,:,i), 'Color',edgecolors(i,:), 'LineWidth',1);
end

xlim([-1*before,after+laser_duration])
ylim([0 100])
xlabel('Time (s)')
ylabel('Percentage (%)')
set(gca,'YTick',0:20:100,'LineWidth',1,'FontSize',12)
[~,lh] = legend([legend_shading(fliplr(plot_order)),legend_mean(fliplr(plot_order)),...
    legend_top(fliplr(plot_order)),legend_bottom(fliplr(plot_order))],...
    [fliplr(state_labels),fliplr(state_labels),...
    fliplr(state_labels),fliplr(state_labels)],'FontSize',13);
legend('boxoff')
box off

lh(13).FaceAlpha = bar_alphas(plot_order(3));
lh(14).FaceAlpha = bar_alphas(plot_order(2));
lh(15).FaceAlpha = bar_alphas(plot_order(1));
lh(13).EdgeAlpha = 0;
lh(14).EdgeAlpha = 0;
lh(15).EdgeAlpha = 0;
lh(16).YData = lh(16).YData + (mean(lh(13).YData) - mean(lh(16).YData));
lh(18).YData = lh(18).YData + (mean(lh(14).YData) - mean(lh(18).YData));
lh(20).YData = lh(20).YData + (mean(lh(15).YData) - mean(lh(20).YData));
lh(22).YData = lh(13).YData(2:3);
lh(24).YData = lh(14).YData(2:3);
lh(26).YData = lh(15).YData(2:3);
lh(28).YData = lh(13).YData([1,4]);
lh(30).YData = lh(14).YData([1,4]);
lh(32).YData = lh(15).YData([1,4]);
for i = 4:12
    lh(i).delete;
end
%% from Franz's code
P = zeros(1,3);
for i=1:3
    if mean(diffs(:,i)) >= 0
        p = length(find(diffs(:,i)>0)) / iters;
        sig = 1 - p;
        if sig == 0, sig = 1/iters; end
    else
        p = length(find(diffs(:,i)<0)) / iters;
        sig = 1 - p;
        if sig == 0, sig = 1/iters; end
    end
    P(i) = sig;
end

disp('p values (pre-laser vs. laser on):')
disp(['REM: ',num2str(P(1))])
disp(['Wake: ',num2str(P(2))])
disp(['NREM: ',num2str(P(3))])

%% transition analysis - preprocessing
binwidth = 60; % width of transition bins in sec
state_bin_new = 20; % resample so each brain state "epoch" is this long, in sec

chunklen = state_bin_new / state_bin;
tp = cell(1,nMice); %zeke's transition probability

allMiceDataBinned = cell(1,nMice); % data binned at new chunk length

for m = 1:nMice % for each mouse
    % resample the time-locked data
    % preallocate
    allMiceDataBinned{m} = zeros(size(allMiceData{m},1), size(allMiceData{m},2)/chunklen);
    
    % downsample the brain state data
    for i = 1:size(allMiceData{m},1) % for each trial
        trialdata = allMiceData{m}(i,:);
        chunk2 = zeros(1,length(trialdata) / chunklen);
        for j = 1:length(chunk2) % for each of the new, larger bins
            chunk = trialdata(((j-1)*chunklen + 1) : (j*chunklen));
            if mode(chunk) == abs(mode(-1*chunk))
                chunk2(j) = mode(chunk);
            else
                if rand > 0.5
                    chunk2(j) = mode(chunk);
                else
                    chunk2(j) = abs(mode(-1*chunk));
                end
            end
        end
        allMiceDataBinned{m}(i,:) = chunk2;
    end
end

%% calculate transitions using binned data
for m = 1:nMice % for each mouse
    tpm = []; % initialize transition probability matrix
    for k = 1:size(allMiceDataBinned{m}, 1) % for each trial
        tpm(:,:,k) = AS_getTransitions(allMiceDataBinned{m}(k,:), state_bin_new, binwidth);
    end
    tp{m} = tpm;
end
tp_mtx = [];
for i = 1:nMice % find avg for each mouse
    tp_mtx(:,:,i) = mean(tp{i},3);
end

stateRows = [1, 2, 3; 4, 5, 6; 7, 8, 9]; % which rows leave which state

tp_avg = mean(tp_mtx,3);
tp_95ci = [];
for j = 1:size(tp_mtx,2) % each timepoint
    for k = 1:3 % each brain state
        allentries = sum(tp_avg(stateRows(k,:),j))+0.0000001; % prevents div/0
        tp_avg(stateRows(k,:),j) = tp_avg(stateRows(k,:),j)/allentries;
    end
end

%% bootstrap
totalpts = (before+laser_duration+after)/binwidth;
% find baseline
baselines = zeros(9,iters);
diffs = zeros(9,iters);
bootstrap_tps = zeros(9,totalpts,iters);

prelsrpts = before/binwidth;
lsrpts = laser_duration / binwidth;

stateRows = [1, 2, 3; 4, 5, 6; 7, 8, 9]; % which rows leave which state

mousetrials = zeros(1,nMice); % number of trials for each mouse
for m = 1:nMice
    mousetrials(m) = size(allMiceDataBinned{m}, 1);
end

wb = waitbar(0,'Bootstrapping...');
for itr = 1:iters
    waitbar(itr/iters,wb);
    tpb = cell(1, nMice); % for bootstrap
    for m = 1:nMice
        j = 1;
        tpm = zeros(9,totalpts,mousetrials(m));
        for k = randi(mousetrials(m),1,mousetrials(m))
            tpm(:,:,j) = AS_getTransitions(allMiceDataBinned{m}(k,:), state_bin_new, binwidth);
            j = j+1;
        end
        tpb{m} = tpm;
    end
    tp_mtx = zeros(9,totalpts,nMice);
    for i = 1:nMice % find avg for each mouse
        tp_mtx(:,:,i) = mean(tpb{i},3);
    end
    
    tp_avgb = mean(tp_mtx,3);
    
    for j = 1:totalpts % each timepoint
        for k = 1:3 % each brain state
            allentries = sum(tp_avgb(stateRows(k,:),j))+0.0000001; % prevents div/0
            tp_avgb(stateRows(k,:),j) = tp_avgb(stateRows(k,:),j)/allentries;
        end
    end
    
    bootstrap_tps(:,:,itr) = tp_avgb;
    baselines(:,itr) = mean(tp_avgb(:,1:prelsrpts),2);
    diffs(:,itr) = mean(tp_avgb(:,(prelsrpts+1):(prelsrpts + lsrpts)),2) - baselines(:,itr);
end
delete(wb);

for i = 1:size(bootstrap_tps,1)
    for j = 1:size(bootstrap_tps,2)
        tp_95ci(i,j,1)=prctile(bootstrap_tps(i,j,:),5);
        tp_95ci(i,j,2)=prctile(bootstrap_tps(i,j,:),95);
    end
end

% from Franz's code
P = zeros(1,9);
for i=1:9
    if mean(diffs(i,:)) >= 0
        p = length(find(diffs(i,:)>0)) / iters;
        sig = 1 - p;
        if sig == 0, sig = 1/iters; end
    else
        p = length(find(diffs(i,:)<0)) / iters;
        sig = 1 - p;
        if sig == 0, sig = 1/iters; end
    end
    P(i) = sig;
end

%% plot transition probabilities
bl_mean = mean(baselines,2);

figure('Color','w','Position',[900,296,811,653])
labels = {'REM \rightarrow REM','REM \rightarrow Wake','REM \rightarrow NREM',...
    'Wake \rightarrow REM','Wake \rightarrow Wake','Wake \rightarrow NREM',...
    'NREM \rightarrow REM','NREM \rightarrow Wake','NREM \rightarrow NREM'};
labelstxt = {'REM -> REM','REM -> Wake','REM -> NREM','Wake -> REM','Wake -> Wake',...
    'Wake -> NREM','NREM -> REM','NREM -> Wake','NREM -> NREM'};
xpos = [.08 .4 .71 .08 .4 .71 .08 .4 .71];
ypos = [.71 .71 .71 .39 .39 .39 .07 .07 .07];

for i = 1:9
    axes('Position',[xpos(i) ypos(i) .25 .25])
    fill([before before+laser_duration before+ laser_duration before],...
        [0 0 1 1],[0 0 1], 'EdgeColor', 'none','FaceAlpha', 0.2)
    hold on, bar(binwidth*(1:size(tp_avg,2))-binwidth/2,tp_avg(i,:),'FaceColor',[.98 .98 .98]);
    hold on, plot([binwidth*(1:size(tp_avg,2))-binwidth/2; binwidth*(1:size(tp_avg,2))-binwidth/2],...
        [tp_95ci(i,:,1);tp_95ci(i,:,2)],'k')
    
    title(labels{i})
    upperlim = max(tp_95ci(i,:,2));
    lowerlim = min(tp_95ci(i,:,1));
    
    if upperlim < .1
        ylim([0 .1])
    else
        if upperlim < .25
            ylim([0 .25])
        else
            if upperlim < .5
                ylim([0 .5])
            else
                if lowerlim > 0.6
                    ylim([floor(lowerlim*10)/10 1])
                else
                    ylim([0 1])
                end
            end
        end
    end
    xlim([0, before+after+laser_duration])
    
    if i == 8
        xlabel('Time (s)')
    end
    if i == 4
        ylabel('Transition probability')
    end
    
    % plot baseline
    hold on, plot(xlim, ones(1,2)*bl_mean(i),'r--')
    box off
    
    if laser_duration == 60
        set(gca,'XTick',0:60:(laser_duration + before+after),'XTickLabel',...
            (0:60:(laser_duration + before+after))-before)
    else
        set(gca,'XTick',[0 120 240 360 480 600],'XTickLabel',[-240 -120 0 120 240 360])
    end
end

% display p values
disp('p values (pre-laser vs. laser on):')
for i = 1:9
    disp([labelstxt{i},' : ',num2str(P(i))])
end
