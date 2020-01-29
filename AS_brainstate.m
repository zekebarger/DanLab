% zeke barger 012920
% plots brain state relative to laser onset
% averages across animals and uses bootstrapping
% requires a mouselist (list of dirlist file locations)

%TODO
% eyfp comparison
% arbitrary time axis

function [] = AS_brainstate()
%% set parameters
% general parameters:
iters = 100; % iterations for generating bootstrap confidence intervals (default 10000)
before = 240; % seconds before laser onset to display (default 240)
after = 240; % seconds after to display (default 240)
laser_duration = 120; % laser duration in sec
state_bin = 2.5; % epoch length for brain state classification
% parameters for plotting brain state over time:
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
% parameters for transition analysis:
binwidth = 60; % width of transition bins in sec
state_bin_new = 20; % resample so each brain state "epoch" is this long, in sec
sig_thresh = 0.05; % significance threshold
% because we downsample the brain states to a lower time resolution, there
% are usually some 'impossible' transitions between R-N and W-R. Setting
% the flags below to 1 will convert these transitions to R-W and N-R,
% respectively.
rn_2_rw = 1;
wr_2_nr = 1;


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
warning('off','MATLAB:load:variableNotFound');
temp = load([fpath,fname],'mouselist'); % load list of dirlists
warning('on','MATLAB:load:variableNotFound');
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
box off

% legend ;)
[~,lh] = legend([legend_shading(fliplr(plot_order)),legend_mean(fliplr(plot_order)),...
    legend_top(fliplr(plot_order)),legend_bottom(fliplr(plot_order))],...
    [fliplr(state_labels),fliplr(state_labels),...
    fliplr(state_labels),fliplr(state_labels)],'FontSize',13);
legend('boxoff')
lh(13).FaceAlpha = bar_alphas(plot_order(3)); % fix shading alpha
lh(14).FaceAlpha = bar_alphas(plot_order(2));
lh(15).FaceAlpha = bar_alphas(plot_order(1));
lh(13).EdgeAlpha = 0; % fix edge alpha
lh(14).EdgeAlpha = 0;
lh(15).EdgeAlpha = 0;
% align lines with shading
lh(16).YData = lh(16).YData + (mean(lh(13).YData) - mean(lh(16).YData));
lh(18).YData = lh(18).YData + (mean(lh(14).YData) - mean(lh(18).YData));
lh(20).YData = lh(20).YData + (mean(lh(15).YData) - mean(lh(20).YData));
lh(22).YData = lh(13).YData(2:3);
lh(24).YData = lh(14).YData(2:3);
lh(26).YData = lh(15).YData(2:3);
lh(28).YData = lh(13).YData([1,4]);
lh(30).YData = lh(14).YData([1,4]);
lh(32).YData = lh(15).YData([1,4]);
for i = 4:12 % remove entry labels we don't need
    lh(i).String='';
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
    if wr_2_nr
        tpm(7,:,:) = tpm(7,:,:) + tpm(4,:,:);
        tpm(4,:,:) = 0;
    end
    if rn_2_rw
        tpm(2,:,:) = tpm(2,:,:) + tpm(3,:,:);
        tpm(3,:,:) = 0;
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
        if wr_2_nr
            tpm(7,:,:) = tpm(7,:,:) + tpm(4,:,:);
            tpm(4,:,:) = 0;
        end
        if rn_2_rw
            tpm(2,:,:) = tpm(2,:,:) + tpm(3,:,:);
            tpm(3,:,:) = 0;
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
labels = {'R \rightarrow R','R \rightarrow W','R \rightarrow N',...
    'W \rightarrow R','W \rightarrow W','W \rightarrow N',...
    'N \rightarrow R','N \rightarrow W','N \rightarrow N'};
labelstxt = {'REM -> REM','REM -> Wake','REM -> NREM','Wake -> REM','Wake -> Wake',...
    'Wake -> NREM','NREM -> REM','NREM -> Wake','NREM -> NREM'};
xpos = [.08 .4 .71 .08 .4 .71 .08 .4 .71];
ypos = [.72 .72 .72 .40 .40 .40 .08 .08 .08];
th = {};
titleFS = 14;
asteriskFS = 35;
tplot_order = [5 4 6 3 1 2 8 7 9];
ac = {};
for i = 1:9
    if (rn_2_rw && i == 3) || (wr_2_nr && i == 4)
        continue
    end
    axes('Position',[xpos(tplot_order(i)) ypos(tplot_order(i)) .24 .20])
    fill([before before+laser_duration before+ laser_duration before],...
        [0 0 1 1],[0 0 1], 'EdgeColor', 'none','FaceAlpha', 0.2)
    hold on, bar(binwidth*(1:size(tp_avg,2))-binwidth/2,tp_avg(i,:),'FaceColor',[.98 .98 .98]);
    hold on, plot([binwidth*(1:size(tp_avg,2))-binwidth/2; binwidth*(1:size(tp_avg,2))-binwidth/2],...
        [tp_95ci(i,:,1);tp_95ci(i,:,2)],'k')
    
    upperlim = max(tp_95ci(i,:,2));
    lowerlim = min(tp_95ci(i,:,1));
    
    if upperlim < .3
        ylim([0 .3])
    else
        ylim([0 1])
    end
    xlim([0, before+after+laser_duration])
    
    if i == 8
        xlabel('Time (s)')
    end
    if i == 5 || (i == 2 || i ==8)
        ylabel('Probability')
    end
    
    % plot baseline
    hold on, plot(xlim, ones(1,2)*bl_mean(i),'r--')
    box off
    
    if laser_duration == 60
        set(gca,'XTick',0:60:(laser_duration + before+after),'XTickLabel',...
            (0:60:(laser_duration + before+after))-before)
    else
        set(gca,'XTick',[0 120 240 360 480 600],...
            'XTickLabel',[-240 -120 0 120 240 360],...
            'FontSize',12)
    end
    
    % title
    title({labels{i};'*'},'Color','w','FontSize',titleFS); % placeholder
    th{i}=get(gca,'Title');
    if mean(diffs(i,:)) > 0 && P(i) < sig_thresh
        ac{i} = 'magenta';
        as = '*';
    else
        if P(i) < sig_thresh
            ac{i} = 'cyan';
            as = '*';
        else
            ac{i} = 'white';
            as = '';
        end
    end
    text(th{i}.Position(1),th{i}.Position(2)*.94,{labels{i};''},...
        'HorizontalAlignment','center','FontSize',titleFS,'Margin',3,'VerticalAlignment','bottom')
    text(th{i}.Position(1),th{i}.Position(2)*.78,{['\fontsize{',num2str(asteriskFS),...
        '}\color{',ac{i},'}',as]},...
        'HorizontalAlignment','center','FontSize',13,'Margin',3,'VerticalAlignment','bottom')
end

% draw a diagram!
axes('Position',[.66 .45 .28 .55])
% draw letters
text(.46, .69, 'N','FontSize',21)
text(.20, .28, 'W','FontSize',21)
text(.72, .28, 'R','FontSize',21)
for i = 1:9 % convert white colors to black
    if strcmp(ac{i},'white')
        ac{i} = 'k';
    end
end
% draw straight arrows
annotation('arrow',[.821 .870],[.758 .680],'HeadStyle','plain','LineWidth',2,...
    'HeadLength',7,'HeadWidth',7,'Color',ac{7});
annotation('arrow',[.847 .760],[.652 .652],'HeadStyle','plain','LineWidth',2,...
    'HeadLength',7,'HeadWidth',7,'Color',ac{2});
annotation('arrow',[.732 .779],[.691 .766],'HeadStyle','plain','LineWidth',2,...
    'HeadLength',7,'HeadWidth',7,'Color',ac{6});
annotation('arrow',[.791 .745],[.753 .680],'HeadStyle','plain','LineWidth',2,...
    'HeadLength',7,'HeadWidth',7,'Color',ac{8});
% draw curved arrows
hold on, circular_arrow(.08, [.50 .82], 90, 280, 1, ac{9}, 9, .08);
hold on, circular_arrow(.08, [.84 .19], 315, 280, 1, ac{1}, 9, .15);
hold on, circular_arrow(.08, [.175 .19], 235, 280, 1, ac{5}, 9, .14);
axis equal
axis off

% legend
hold on, p1=plot([10 10],[10 10],'m','LineWidth',2);
hold on, p2=plot([10 10],[10 10],'c','LineWidth',2);
hold on, p3=plot([10 10],[10 10],'k','LineWidth',2);
legend([p1 p2 p3],{'Significant increase','Significant decrease','No change'},...
    'Position',[0.6775,0.4589,0.2416,0.1041],'FontSize',12)
legend('boxoff')

xlim([0 1])
ylim([0 1])

% display p values
disp('p values (pre-laser vs. laser on):')
for i = 1:9
    disp([labelstxt{i},' : ',num2str(P(i))])
end


function circular_arrow(radius, centre, arrow_angle, angle, direction, colour, head_size, qq)
% by Zac Giles https://www.mathworks.com/matlabcentral/fileexchange/59917-circular_arrow
arrow_angle = deg2rad(arrow_angle); % Convert angle to rad
angle = deg2rad(angle); % Convert angle to rad
xc = centre(1);
yc = centre(2);
x_temp = centre(1) + radius;
y_temp = centre(2);
x1 = (x_temp-xc)*cos(arrow_angle+angle/2) - ...
        (y_temp-yc)*sin(arrow_angle+angle/2) + xc;
x2 = (x_temp-xc)*cos(arrow_angle-angle/2) - ...
        (y_temp-yc)*sin(arrow_angle-angle/2) + xc;
x0 = (x_temp-xc)*cos(arrow_angle) - ...
        (y_temp-yc)*sin(arrow_angle) + xc;
y1 = (x_temp-xc)*sin(arrow_angle+angle/2) + ...
        (y_temp-yc)*cos(arrow_angle+angle/2) + yc;
y2 = (x_temp-xc)*sin(arrow_angle-angle/2) + ... 
        (y_temp-yc)*cos(arrow_angle-angle/2) + yc;
y0 = (x_temp-xc)*sin(arrow_angle) + ... 
        (y_temp-yc)*cos(arrow_angle) + yc;
i = 1;
P1 = struct([]);
P2 = struct([]);
P1{1} = [x1;y1]; % Point 1 - 1
P1{2} = [x2;y2]; % Point 1 - 2
P2{1} = [x0;y0]; % Point 2 - 1
P2{2} = [x0;y0]; % Point 2 - 1
centre = [xc;yc]; % guarenteeing centre is the right dimension
n = 1000; % The number of points in the arc
v = struct([]);
while i < 3
    v1 = P1{i}-centre;
    v2 = P2{i}-centre;
    c = det([v1,v2]); % "cross product" of v1 and v2
    a = linspace(0,atan2(abs(c),dot(v1,v2)),n); % Angle range
    v3 = [0,-c;c,0]*v1; % v3 lies in plane of v1 and v2 and is orthog. to v1
    v{i} = v1*cos(a)+((norm(v1)/norm(v3))*v3)*sin(a); % Arc, center at (0,0)
    plot(v{i}(1,:)+xc,v{i}(2,:)+yc,'Color', colour,'LineWidth',2) % Plot arc, centered at P0
    i = i + 1;
end
position = struct([]);
if direction == 1
%     position{1} = [x2 y2 x2-(v{2}(1,2)+xc) y2-(v{2}(2,2)+yc)];
    position{1} = [(v{2}(1,30)+xc), (v{2}(2,30)+yc),...
        v{2}(1,10)-v{2}(1,n*qq), v{2}(2,10)-v{2}(2,n*qq)];
elseif direction == -1
    position{1} = [x1 y1 x1-(v{1}(1,2)+xc) y1-(v{1}(2,2)+yc)];
elseif direction == 2
    position{1} = [x2 y2 x2-(v{2}(1,2)+xc) y2-(v{2}(2,2)+yc)];
    position{2} = [x1 y1 x1-(v{1}(1,2)+xc) y1-(v{1}(2,2)+yc)];  
elseif direction == 0
    % Do nothing
else
    error('direction flag not 1, -1, 2 or 0.');
end
i = 1;
while i < abs(direction) + 1
    h=annotation('arrow'); % arrow head
    set(h,'parent', gca, 'position', position{i}, ...
        'HeadLength', head_size, 'HeadWidth', head_size,...
        'HeadStyle', 'plain', 'linestyle','none','Color', colour);
    i = i + 1;
end
