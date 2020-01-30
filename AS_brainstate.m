% zeke barger 013020
% plots brain state relative to laser onset
% averages across animals and uses bootstrapping
% requires a mouselist (list of dirlist file locations) or dirlist

function [] = AS_brainstate()
%% set parameters
% general parameters:
iters = 10000; % iterations for generating bootstrap confidence intervals (default 10000)
before = 240; % seconds before laser onset to display (default 240)
after = 240; % seconds after to display (default 240)
laser_duration = 120; % laser duration in sec
epoch_length = 2.5; % epoch length for brain state classification
% parameters for plotting brain state over time:
% order to plot the lines (last entry is plotted on top). 1=R, 2=W, 3=N
plot_order = [1 2 3];
% line colors (REM, wake, NREM)
colors = [43 160 220;
    129 131 132
    255 194 61] ./ 255;
% error bar colors
error_colors = [43 160 220;
    81 82 82
    255 194 61] ./ 255;
% color for edge of error bars
edge_colors = [150 207 238;
    194 195 196;
    244 212 153] ./ 255;
% transparency of error bars
bar_alphas = [.4 .4 .4];
% parameters for transition analysis:
bin_width = 60; % width of transition bins in sec
transition_epoch_length = 20; % resample so each brain state "epoch" is this long, in sec
sig_thresh = 0.05; % significance threshold
% because we downsample the brain states to a lower time resolution, there
% are usually some 'impossible' transitions between R-N and W-R. Setting
% the flags below to 1 will convert these transitions to R-W and N-R,
% respectively. If both are not set to 1, the transition diagram will not
% be shown.
rn_2_rw = 1;
wr_2_nr = 1;


% no need to edit these lines
% epochs to plot before, during and after laser
before_epochs = round(before/epoch_length); % number of bins before laser
after_epochs = round(after/epoch_length); % number of bins after laser
laser_epochs = round(laser_duration / epoch_length); % number of bins during laser

fig_positions = {[137   562   560   420],...
    [703   328   811   653];...
    [137    55   560   420],...
    [703     7   811   653]};
%% load and process data
% have user select the mouselist
[fname,fpath,~] = uigetfile('*.mat','Select mouselist for TREATMENT condition, or a dirlist');
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
    treatment_list = temp.dirlist;
    dirlist_mode = 1;
else
    treatment_list = temp.mouselist;
    dirlist_mode = 0;
end

if isempty(treatment_list) % nothing in the list
    disp('the list is empty');
    return
end

% then, if a mouselist was selected, provide the option to choose a CONTROL
% list for comparison
control_list = {};
if ~dirlist_mode
    [fname,fpath,~] = uigetfile([fpath,'*.mat'],'Select mouselist for CONTROL condition, or CANCEL');
    % check if something was selected
    if ischar(fname) && ~isempty(whos('-file',[fpath,fname], 'mouselist'))
        temp = load([fpath,fname],'mouselist'); % load list of dirlists
        control_list = temp.mouselist;
    else
        disp('No mouselist was selected.');
    end
end

% put treatment and control mouselists into one array
all_mouselists = {treatment_list,control_list};
% determine whether to analyze both mouselists
num_lists = 1 + ~isempty(control_list);

% process each list
experiment_types = {'TREATMENT','CONTROL'};
for this_list = 1:num_lists
    mouselist = all_mouselists{this_list};
    % load data
    num_mice = length(mouselist); % how many mice there are
    labels_all_mice = cell(1, num_mice); % laser locked brain states, one cell per mouse
    filenames_all_mice = cell(1, num_mice); % locations of files
    num_days_all_mice = zeros(1,num_mice); % how many days there are for each mouse
    % if we are dealing with a mouselist
    if ~dirlist_mode
        % check that the contents of the mouselist are ok, and gather filenames
        for m = 1:num_mice % for each mouse
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
            filenames_all_mice{m} = fileNames;
            num_days_all_mice(m) = length(dirlist);
        end
    else % we are dealing with a dirlist
        [problemString, fileNames] = AS_checkList(mouselist, {'labels','laser'});
        % if something is wrong with this list
        if ~isempty(problemString)
            % display the problem and quit
            disp(problemString)
            return
        end
        for m = 1:num_mice % for each recording
            filenames_all_mice{m} = {fileNames{m}};
            num_days_all_mice(m) = 1; % just one "day"
        end
    end
    
    wb = waitbar(0,'Loading...');
    for m = 1:num_mice % for each mouse
        waitbar(m/num_mice,wb); % update progress bar
        labels_all_mice{m} = []; % initialize array of trial info
        for d = 1:num_days_all_mice(m) % for each day
            % get laser-locked brain states
            labels_all_mice{m} = [labels_all_mice{m}; AS_getLaserStates(filenames_all_mice{m}{d}{1},...
                filenames_all_mice{m}{d}{2},before_epochs,laser_epochs,after_epochs, epoch_length)];
        end
    end
    delete(wb)
    
    % for each mouse, average across days
    mouseAvgs = []; % state, time, mouse
    for i = 1:3 % state
        for j = 1:num_mice % mouse
            mouseAvgs(i,:,j) = mean(labels_all_mice{j}==i);
        end
    end
    
    %% find mean and CIs
    trialsPerMouse = zeros(1,num_mice); % number of trials for each mouse
    for i = 1:num_mice
        trialsPerMouse(i) = size(labels_all_mice{i},1);
    end
    % totaltrials = sum(mousetrials); % total number of trials in the whole experiment
    trialdur = size(labels_all_mice{1},2); % number of time bins in each trial
    all_mice_boot_avg = zeros(iters,trialdur,3); % iter by time by state bootstrapped likelihood
    
    wb = waitbar(0,'Bootstrapping...');
    state_diffs{this_list} = zeros(iters, 3); % pre vs laser difference
    pre_lsr_pts = before / epoch_length;
    lsr_pts = laser_duration / epoch_length;
    for itr = 1:iters
        waitbar(itr/iters,wb);
        mouse_boot_avg = zeros(num_mice,trialdur,3);
        for m = 1:num_mice % get data randomly from each mouse
            bsdata = labels_all_mice{m}(randi(trialsPerMouse(m),1,trialsPerMouse(m)),:);
            % average across trials for this mouse
            for k = 1:3 % state
                mouse_boot_avg(m,:,k) = mean(bsdata==k);
            end
        end
        
        all_mice_boot_avg(itr,:,:) = mean(mouse_boot_avg,1);
        state_diffs{this_list}(itr,:) = squeeze(mean(all_mice_boot_avg(itr,1:pre_lsr_pts,:))) -...
            squeeze(mean(all_mice_boot_avg(itr,(1:lsr_pts)+pre_lsr_pts,:)));
        
    end
    delete(wb);
    
    CIs = prctile(all_mice_boot_avg,[2.5 97.5],1);
    CIs = CIs*100;
    grandMean = 100*mean(mouseAvgs,3);
    
    %% plot
    state_labels = {'REM','Wake','NREM'};
    state_labels = state_labels(plot_order);
    
    figure('Color','w','Position',fig_positions{this_list,1})
    fill([0 laser_duration laser_duration 0],[0 0 100 100],[154 154 255]./255,...
        'EdgeAlpha',0,'FaceAlpha',.5);
    
    t = (1:(after_epochs+before_epochs+laser_epochs))*epoch_length -...
        before_epochs*epoch_length - epoch_length/2;
    t_flip = cat(2,t,fliplr(t));
    
    legend_mean = [];
    legend_shading = [];
    legend_top = [];
    legend_bottom = [];
    for i = plot_order
        hold on, legend_shading(i) = fill(t_flip, cat(2,CIs(1,:,i),fliplr(CIs(2,:,i))),...
            error_colors(i,:), 'EdgeAlpha', 0,'FaceAlpha', bar_alphas(i),'LineWidth',1.5);
        hold on, legend_mean(i) = plot(t, grandMean(i,:), 'Color',colors(i,:), 'LineWidth',1.5);
        hold on, legend_top(i) = plot(t, CIs(2,:,i), 'Color',edge_colors(i,:), 'LineWidth',1);
        hold on, legend_bottom(i) = plot(t, CIs(1,:,i), 'Color',edge_colors(i,:), 'LineWidth',1);
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
    % from Franz's code: bootstrap significance testing
    P = zeros(1,3);
    for i=1:3
        if mean(state_diffs{this_list}(:,i)) >= 0
            sig = 1 - sum(state_diffs{this_list}(:,i)>0) / iters;
            if sig == 0, sig = 1/iters; end
        else
            sig = 1 - sum(state_diffs{this_list}(:,i)<0) / iters;
            if sig == 0, sig = 1/iters; end
        end
        P(i) = sig;
    end
    
    disp([experiment_types{this_list},' p values (pre-laser vs. laser on):'])
    disp(['REM: ',num2str(P(1))])
    disp(['Wake: ',num2str(P(2))])
    disp(['NREM: ',num2str(P(3))])
    
    %% transition analysis - preprocessing
    % find out how much to downsample the brain state labels
    downsample_epochs = transition_epoch_length / epoch_length;
    all_tp = cell(1,num_mice); % transition probability
    all_labels_binned = cell(1,num_mice); % data binned at new epoch length
    
    for m = 1:num_mice % for each mouse
        % resample the time-locked data
        % preallocate
        all_labels_binned{m} = zeros(size(labels_all_mice{m},1),...
            size(labels_all_mice{m},2)/downsample_epochs);
        
        % downsample the brain state data
        for i = 1:size(labels_all_mice{m},1) % for each trial
            trialdata = labels_all_mice{m}(i,:);
            new_trial = zeros(1,length(trialdata) / downsample_epochs);
            for j = 1:length(new_trial) % for each of the new, larger epochs
                segment = trialdata(((j-1)*downsample_epochs + 1) : (j*downsample_epochs));
                if mode(segment) == abs(mode(-1*segment))
                    new_trial(j) = mode(segment);
                else
                    if rand > 0.5
                        new_trial(j) = mode(segment);
                    else
                        new_trial(j) = abs(mode(-1*segment));
                    end
                end
            end
            all_labels_binned{m}(i,:) = new_trial;
        end
    end
    
    %% calculate transitions using binned data
    for m = 1:num_mice % for each mouse
        mouse_tp = []; % initialize transition probability matrix
        for k = 1:size(all_labels_binned{m}, 1) % for each trial
            mouse_tp(:,:,k) = AS_getTransitions(all_labels_binned{m}(k,:),...
                transition_epoch_length, bin_width);
        end
        if wr_2_nr
            mouse_tp(7,:,:) = mouse_tp(7,:,:) + mouse_tp(4,:,:);
            mouse_tp(4,:,:) = 0;
        end
        if rn_2_rw
            mouse_tp(2,:,:) = mouse_tp(2,:,:) + mouse_tp(3,:,:);
            mouse_tp(3,:,:) = 0;
        end
        all_tp{m} = mouse_tp;
    end
    tp_all_means = [];
    for i = 1:num_mice % find avg for each mouse
        tp_all_means(:,:,i) = mean(all_tp{i},3);
    end
    
    stateRows = [1, 2, 3; 4, 5, 6; 7, 8, 9]; % which rows leave which state
    
    tp_grand_mean = mean(tp_all_means,3);
    tp_95ci = [];
    for j = 1:size(tp_all_means,2) % each timepoint
        for k = 1:3 % each brain state
            allentries = sum(tp_grand_mean(stateRows(k,:),j))+0.0000001; % prevents div/0
            tp_grand_mean(stateRows(k,:),j) = tp_grand_mean(stateRows(k,:),j)/allentries;
        end
    end
    
    %% bootstrap
    totalpts = (before+laser_duration+after)/bin_width;
    % find baseline
    baselines = zeros(9,iters);
    transition_diffs{this_list} = zeros(9,iters);
    bootstrap_tps = zeros(9,totalpts,iters);
    
    prelsrpts = before/bin_width;
    lsrpts = laser_duration / bin_width;
    
    stateRows = [1, 2, 3; 4, 5, 6; 7, 8, 9]; % which rows leave which state
    
    mousetrials = zeros(1,num_mice); % number of trials for each mouse
    for m = 1:num_mice
        mousetrials(m) = size(all_labels_binned{m}, 1);
    end
    
    wb = waitbar(0,'Bootstrapping...');
    for itr = 1:iters
        waitbar(itr/iters,wb);
        tp_boot = cell(1, num_mice); % for bootstrap
        for m = 1:num_mice
            j = 1;
            mouse_tp = zeros(9,totalpts,mousetrials(m));
            for k = randi(mousetrials(m),1,mousetrials(m))
                mouse_tp(:,:,j) = AS_getTransitions(all_labels_binned{m}(k,:),...
                    transition_epoch_length, bin_width);
                j = j+1;
            end
            if wr_2_nr
                mouse_tp(7,:,:) = mouse_tp(7,:,:) + mouse_tp(4,:,:);
                mouse_tp(4,:,:) = 0;
            end
            if rn_2_rw
                mouse_tp(2,:,:) = mouse_tp(2,:,:) + mouse_tp(3,:,:);
                mouse_tp(3,:,:) = 0;
            end
            tp_boot{m} = mouse_tp;
        end
        tp_all_means = zeros(9,totalpts,num_mice);
        for i = 1:num_mice % find avg for each mouse
            tp_all_means(:,:,i) = mean(tp_boot{i},3);
        end
        
        tp_avg_boot = mean(tp_all_means,3);
        
        for j = 1:totalpts % each timepoint
            for k = 1:3 % each brain state
                allentries = sum(tp_avg_boot(stateRows(k,:),j))+0.0000001; % prevents div/0
                tp_avg_boot(stateRows(k,:),j) = tp_avg_boot(stateRows(k,:),j)/allentries;
            end
        end
        
        bootstrap_tps(:,:,itr) = tp_avg_boot;
        baselines(:,itr) = mean(tp_avg_boot(:,1:prelsrpts),2);
        transition_diffs{this_list}(:,itr) = mean(tp_avg_boot(:,...
            (prelsrpts+1):(prelsrpts + lsrpts)),2) - baselines(:,itr);
    end
    delete(wb);
    
    for i = 1:size(bootstrap_tps,1)
        for j = 1:size(bootstrap_tps,2)
            tp_95ci(i,j,1)=prctile(bootstrap_tps(i,j,:),5);
            tp_95ci(i,j,2)=prctile(bootstrap_tps(i,j,:),95);
        end
    end
    
    P = zeros(1,9);
    for i=1:9
        if mean(transition_diffs{this_list}(i,:)) >= 0
            sig = 1 - sum(transition_diffs{this_list}(i,:)>0) / iters;
            if sig == 0, sig = 1/iters; end
        else
            sig = 1 - sum(transition_diffs{this_list}(i,:)<0) / iters;
            if sig == 0, sig = 1/iters; end
        end
        P(i) = sig;
    end
    
    %% plot transition probabilities
    bl_mean = mean(baselines,2);
    
    figure('Color','w','Position',fig_positions{this_list,2})
    labels = {'R \rightarrow R','R \rightarrow W','R \rightarrow N',...
        'W \rightarrow R','W \rightarrow W','W \rightarrow N',...
        'N \rightarrow R','N \rightarrow W','N \rightarrow N'};
    labelstxt = {'REM -> REM','REM -> Wake','REM -> NREM','Wake -> REM','Wake -> Wake',...
        'Wake -> NREM','NREM -> REM','NREM -> Wake','NREM -> NREM'};
    xpos = [.08 .4 .71 .08 .4 .71 .08 .4 .71];
    ypos = [.72 .72 .72 .40 .40 .40 .08 .08 .08];
    th = {}; % title handles
    titleFS = 14; % title font size
    asteriskFS = 35; % asterisk font size
    tplot_order = [5 4 6 3 1 2 8 7 9]; % transition plot order
    ac = {}; % asterisk color
    for i = 1:9
        if (rn_2_rw && i == 3) || (wr_2_nr && i == 4)
            continue
        end
        axes('Position',[xpos(tplot_order(i)) ypos(tplot_order(i)) .24 .20])
        fill([before before+laser_duration before+ laser_duration before],...
            [0 0 1 1],[0 0 1], 'EdgeColor', 'none','FaceAlpha', 0.2)
        hold on, bar(bin_width*(1:size(tp_grand_mean,2))-bin_width/2,tp_grand_mean(i,:),...
            'FaceColor',[.98 .98 .98]);
        hold on, plot([bin_width*(1:size(tp_grand_mean,2))-bin_width/2; ...
            bin_width*(1:size(tp_grand_mean,2))-bin_width/2],...
            [tp_95ci(i,:,1);tp_95ci(i,:,2)],'k')
        
        upperlim = max(tp_95ci(i,:,2));
        if upperlim < .3
            ylim([0 .3])
        else
            ylim([0 1])
        end
        xlim([0, before+after+laser_duration])
        
        if i == 8
            xlabel('Time (s)')
        end
        if i == 5 || (i == 2 || i == 8)
            ylabel('Probability')
        end
        
        % plot baseline
        hold on, plot(xlim, ones(1,2)*bl_mean(i),'r--')
        box off
        
        if laser_duration < 120
            tick_interval = 60;
        else
            tick_interval = 120;
        end
        set(gca,'XTick',0:tick_interval:(laser_duration + before+after),'XTickLabel',...
            (0:tick_interval:(laser_duration + before+after))-before,...
            'FontSize',12)

        % title
        title({labels{i};'*'},'Color','w','FontSize',titleFS); % placeholder
        th{i}=get(gca,'Title');
        if mean(transition_diffs{this_list}(i,:)) > 0 && P(i) < sig_thresh
            ac{i} = 'magenta';
            as = '*'; % asterisk string
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
    if rn_2_rw && wr_2_nr
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
        hold on, circular_arrow(.08, [.50 .82], 90, 280, ac{9}, 9, .08);
        hold on, circular_arrow(.08, [.84 .19], 315, 280, ac{1}, 9, .15);
        hold on, circular_arrow(.08, [.175 .19], 235, 280, ac{5}, 9, .14);
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
    end
    
    for i = 1:9
        disp([labelstxt{i},' : ',num2str(P(i))])
    end
end

% compare brain state likelihood between tx and ctrl
if num_lists == 2 % if we had both treatment and control experiments
    diffs = state_diffs{2} - state_diffs{1};
    P = zeros(1,3);
    for i=1:3
        if mean(diffs(:,i)) >= 0
            sig = 1 - sum(diffs(:,i)>0) / iters;
            if sig == 0, sig = 1/iters; end
        else
            sig = 1 - sum(diffs(:,i)<0) / iters;
            if sig == 0, sig = 1/iters; end
        end
        P(i) = sig;
    end
    disp([experiment_types{1},' vs. ',experiment_types{2},':'])
    disp('p values (pre-laser vs. laser on):')
    disp(['REM: ',num2str(P(1))])
    disp(['Wake: ',num2str(P(2))])
    disp(['NREM: ',num2str(P(3))])
end

function circular_arrow(radius, centre, arrow_angle, angle, colour, head_size, qq)
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
%     position{1} = [x2 y2 x2-(v{2}(1,2)+xc) y2-(v{2}(2,2)+yc)];
position{1} = [(v{2}(1,30)+xc), (v{2}(2,30)+yc),...
    v{2}(1,10)-v{2}(1,n*qq), v{2}(2,10)-v{2}(2,n*qq)];

h=annotation('arrow'); % arrow head
set(h,'parent', gca, 'position', position{1}, ...
    'HeadLength', head_size, 'HeadWidth', head_size,...
    'HeadStyle', 'plain', 'linestyle','none','Color', colour);
