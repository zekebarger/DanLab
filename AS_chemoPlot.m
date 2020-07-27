% Zeke Barger, 072720
% plot brain states following chemogenetic manipulation
% NOTE: the saline and CNO mouselists should point to the same mice,
% in the same order!
% Significance is calculated using a two-way repeated measures ANOVA
% with Bonferroni correction followed by paired t-tests
% *: p < 0.05, **: p < 0.01, ***: p < 0.001
% It's assumed that each recording starts right after the injection

% TODO: stats on duration distributions?

function [] = AS_chemoPlot()
%% user-defined parameters - just edit in this section
% how many hrs of data to take from the beginning of each recording
hours_to_plot = 3;
epoch_length = 2.5; % in seconds
mins_per_bin = 20; % minutes of brain states to avg for each plotted timepoint
recording_delay = 5; % minutes between injection and recording start
plot_order = [3 2 1]; % 1=R, 2=W, 3=N. So, to plot wake, then NREM, then REM,
%                       set this to [2 3 1]

% the two parameters below are the times, in minutes, between which to
% collect data for the plots of bout duration & number (a bout is an 
% uninterrupted period of some brain state)and bout numbers
bout_start_time = 0; % time to start, in minutes
bout_end_time = 180; % time to end, in minutes

% this specifies the ranges of bout lengths to count for each animal
% when looking at bout duration
% example: to only track bouts of 0-10sec, 11-30sec, or 31+sec, 
% you would put [0 10 30]
bout_length_ranges = [0 5 15 60 180]; 

% when calculating the delay until the first NREM bout, this is the minimum
% length for that bout in seconds
min_1st_nrem_duration = 20; 

%% basic plotting info
titles = {'REM','Wake','NREM'};
colors = [43 160 220;
    129 131 132
    255 194 61] ./ 255; % rem wake nrem

%% calculate some things for plotting bout information
% convert from hours/minutes to epochs
first_histo_epoch = max([round(bout_start_time*60/epoch_length), 1]); % epoch to start histogram
last_histo_epoch = round(bout_end_time*60/epoch_length); % epoch to stop histogram
min_1st_nrem_epochs = ceil(min_1st_nrem_duration / epoch_length);

%% get brain state data
% number of epochs to take - don't change this
total_epochs = round(hours_to_plot*60*60/epoch_length); % epochs from beginning
epochs_per_bin = mins_per_bin*60/epoch_length; % bins per timepoint

[sFileName,sPathName,~] = uigetfile('*.mat','Choose saline mouselist');
[cFileName,cPathName,~] = uigetfile([sPathName,'*.mat'],'Choose CNO mouselist');

% check if something was selected
if ~ischar(sFileName) || ~ischar(cFileName)
    disp('no list selected');
    return
end

% load mouselists
sList = load([sPathName,sFileName],'mouselist');
cList = load([cPathName,cFileName],'mouselist');

if isempty(sList.mouselist) || isempty(cList.mouselist) % nothing in the list
    disp('empty mouselist');
    return
end

% load brain state data and bouts
[sData, sBouts, sLatency] = loadData(sList.mouselist, total_epochs,...
    epochs_per_bin, first_histo_epoch, last_histo_epoch, min_1st_nrem_epochs);
[cData, cBouts, cLatency] = loadData(cList.mouselist, total_epochs,...
    epochs_per_bin, first_histo_epoch, last_histo_epoch, min_1st_nrem_epochs);
sData = sData*100;
cData = cData*100;

% if there was some problem
if isempty(sData) || isempty(cData)
    return
end

% average brain states across mice
sAvg = mean(sData,3);
cAvg = mean(cData,3);

% find variability
sSEM = std(sData,[],3)/sqrt(size(sData,3));
cSEM = std(cData,[],3)/sqrt(size(cData,3));

%% plot brain state over time
% set up time axis
t = (mins_per_bin * (1:(total_epochs/epochs_per_bin)) + recording_delay - mins_per_bin/2) / 60 ;
figure('Color','w','Position',[458 408 1070 268])
% keep track of how many plots we've made
plot_number = 1;
for i = plot_order %  for each state
    subplot(1,3,plot_number)
    plot_handles_c = []; % hold handles for when we make the legend
    plot_handles_s = [];
    % plot mean and sem
    plot_handles_s(1) = errorbar(t,sAvg(i,:),sSEM(i,:),'o-','Color',colors(i,:),...
        'CapSize',0,'MarkerFaceColor','w','LineWidth',1);
    hold on, plot_handles_c(1) = errorbar(t,cAvg(i,:),cSEM(i,:),'o-','Color',colors(i,:),...
        'MarkerFaceColor',colors(i,:),'LineWidth',1,'CapSize',0);
    
    if i == 1
        % create legend, by first drawing placeholder lines
        hold on, plot_handles_s(2) = errorbar([-1 -2],[0 0],[0 0],'o-','Color',...
            colors(2,:),'CapSize',0,'MarkerFaceColor','w','LineWidth',1);
        hold on, plot_handles_c(2) = errorbar([-1 -2],[0 0],[0 0],'o-','Color',colors(2,:),...
            'MarkerFaceColor',colors(2,:),'LineWidth',1,'CapSize',0);
        hold on, plot_handles_s(3) = errorbar([-1 -2],[0 0],[0 0],'o-','Color',...
            colors(3,:),'CapSize',0,'MarkerFaceColor','w','LineWidth',1);
        hold on, plot_handles_c(3) = errorbar([-1 -2],[0 0],[0 0],'o-','Color',colors(3,:),...
            'MarkerFaceColor',colors(3,:),'LineWidth',1,'CapSize',0);
        [~,lh] = legend([plot_handles_s(plot_order(1)),plot_handles_c(plot_order(1)),...
            plot_handles_s(plot_order(2)),plot_handles_c(plot_order(2)),...
            plot_handles_s(plot_order(3)),plot_handles_c(plot_order(3))],...
            {'','','','','Saline','CNO'},'Location','north','Box','off',...
            'FontSize',12,'NumColumns',3,'Position',[0.8,0.72,0.095,0.552]);
        % get objects in the legend and rearrange them ;)
        all_obj = cell(1,6);
        for ob = 1:6
            all_obj{ob} = lh(ob+6).Children.Children;
        end
        cno_x1 = all_obj{6}(1).XData;
        cno_y1 = all_obj{6}(1).YData;
        cno_x2 = all_obj{6}(2).XData;
        cno_y2 = all_obj{6}(2).YData;
        sal_x1 = all_obj{5}(1).XData;
        sal_y1 = all_obj{5}(1).YData;
        sal_x2 = all_obj{5}(2).XData;
        sal_y2 = all_obj{5}(2).YData;
        x_shift = 1.1*diff(cno_x2);
        all_obj{1}(1).YData = sal_y1;
        all_obj{1}(2).YData = sal_y2;
        all_obj{1}(1).XData = all_obj{1}(1).XData - 2*x_shift;
        all_obj{1}(2).XData = all_obj{1}(2).XData - 2*x_shift;
        all_obj{2}(1).YData = cno_y1;
        all_obj{2}(2).YData = cno_y2;
        all_obj{2}(1).XData = all_obj{2}(1).XData - 2*x_shift;
        all_obj{2}(2).XData = all_obj{2}(2).XData - 2*x_shift;
        all_obj{3}(1).YData = sal_y1;
        all_obj{3}(2).YData = sal_y2;
        all_obj{3}(1).XData = all_obj{3}(1).XData - x_shift;
        all_obj{3}(2).XData = all_obj{3}(2).XData - x_shift;
        all_obj{4}(1).YData = cno_y1;
        all_obj{4}(2).YData = cno_y2;
        all_obj{4}(1).XData = all_obj{4}(1).XData - x_shift;
        all_obj{4}(2).XData = all_obj{4}(2).XData - x_shift;
    end
    
    if plot_number==1 % show y label
        ylabel('Percentage (%)')
    end
    title(titles{i},'Color',colors(i,:),'FontSize',13)
    xlabel('Time (h)')
    xlim([0 hours_to_plot])
    set(gca,'FontSize',12)
    box off
    
    yl = ylim;
    if yl(2) > 20
        yl(2) = 100;
    else
        yl(2) = 20;
    end
    ylim([0 yl(2)])
    % check that number of mice and timepoints are equal
    % (the timepoints should be, for sure. mice # could
    % vary if there's a problem with the data).
    
    if size(sData,2) ~= size(cData,2) || size(sData,3) ~= size(cData,3)
        disp('different number of mice in each condition.')
        disp('cannot calculate statistics')
        continue
    end
    
    % perform two-way repeated measures ANOVA to determine whether there is
    % a significant interaction of treatment and time.
    % first, arrange the data properly
    y = zeros(size(sData,2)*size(sData,3)*2, 1); % each measurement
    subject = y; % subject id
    g1 = y; % treatment
    g2 = y; % time
    ii = 1; % just a counter for how many datapoints we've collected
    for j = 1:size(sData,2) % for each timepoint
        for k = 1:size(sData,3) % for each subject
            y(ii) = sData(i,j,k); % get the measurement
            subject(ii) = k; % store the subject id
            g1(ii) = 1; % set treatment to 1
            g2(ii) = j; % store timepoint
            ii = ii+1;
            % repeat for CNO
            y(ii) = cData(i,j,k); % get the measurement
            subject(ii) = k; % store the subject id
            g1(ii) = 2; % set treatment to 2
            g2(ii) = j; % store timepoint
            ii = ii+1;
        end
    end
    % run the ANOVA
    a = rm_anova2(y, subject, g1, g2, {'treatment','time'});
    df1 = num2str(a{4,3});
    df2 = num2str(a{7,3});
    F = num2str(a{4,5});
    pval = a{4,6};
    pval_string = num2str(pval);
    disp([titles{i},' ANOVA (treatment x time): F(',df1,',',df2,') = ',F,', p = ',pval_string])
    % if there was a significant treatment x time interaction, run post-hoc
    % tests at each timepoint
    if pval <= (0.05 / 3) % Bonferroni correction
        for x = 1:length(t)
            [~,p] = ttest(squeeze(sData(i,x,:)),squeeze(cData(i,x,:)));       
            if p > 0.05 || isnan(p)
                continue
            else
                mark = p_val_to_marker(p);
            end
            yl = ylim;
            y=max([sAvg(i,x)+sSEM(i,x), cAvg(i,x)+cSEM(i,x)])+.028*diff(yl);
            text(t(x), y,mark,'HorizontalAlignment','center','FontSize',14)         
        end
    end 
    plot_number = plot_number + 1;
end

%% bout length
% first, concatenate lists of bout lengths
n_mice = size(sBouts,2); % find number of mice
if n_mice ~= size(cBouts,2) %(should be the same in each condition)
    disp('number of mice should be the same in each condition');
    disp('and they should be in the same order!');
    return
end

sBouts_concatenated = cell(size(sBouts));
cBouts_concatenated = cell(size(cBouts));
for i = 1:3
    for j = 1:n_mice
        sBouts_concatenated{i,j} = [];
        for k = 1:length(sBouts{i,j})
            sBouts_concatenated{i,j} = [sBouts_concatenated{i,j}, sBouts{i,j}{k}];
        end

        cBouts_concatenated{i,j} = [];
        for k = 1:length(cBouts{i,j})
            cBouts_concatenated{i,j} = [cBouts_concatenated{i,j}, cBouts{i,j}{k}];
        end
    end
end
% then figure out how to make the histograms
nBins = length(bout_length_ranges); % get number of bins
% preallocate for the histograms from each mouse
sHistos = zeros(3,nBins,n_mice);
cHistos = zeros(3,nBins,n_mice);
% calculate histogram for each mouse
for i = 1:3 % for each state
    for j = 1:n_mice
        [sHistos(i,:,j),~] = histcounts(sBouts_concatenated{i,j}, ...
            [bout_length_ranges, Inf],'Normalization','probability');
        [cHistos(i,:,j),~] = histcounts(cBouts_concatenated{i,j}, ...
            [bout_length_ranges, Inf],'Normalization','probability');
    end
end

% average across mice
sHistAvg = mean(sHistos,3);
cHistAvg = mean(cHistos,3);

% find variability (SEM)
sHistSEM = std(sHistos,[],3)/sqrt(n_mice);
cHistSEM = std(cHistos,[],3)/sqrt(n_mice);

% then average durations over days, then mice for each state
sDuration = zeros(3,n_mice);
cDuration = zeros(3,n_mice);
for i = 1:3 % for each state
    for j = 1:n_mice % for each mouse
        sMeans = zeros(1,length(sBouts{i,j}));
        cMeans = zeros(1,length(cBouts{i,j}));
        for k = 1:length(sBouts{i,j}) % for each recording
            sMeans(k) = mean(sBouts{i,j}{k});
        end
        for k = 1:length(cBouts{i,j}) % for each recording
            cMeans(k) = mean(cBouts{i,j}{k});
        end
        % average across days
        sDuration(i,j) = mean(sMeans)*epoch_length;
        cDuration(i,j) = mean(cMeans)*epoch_length;
    end
end

% average across mice
sDurationAvg = mean(sDuration,2);
cDurationAvg = mean(cDuration,2);

% find variability (SEM)
sDurationSEM = std(sDuration,[],2)/sqrt(n_mice);
cDurationSEM = std(cDuration,[],2)/sqrt(n_mice);

%% plot simple bout length bar plot
figure('Color','w','Position',[1242 393 347 268])
x = [1 2; 4 5; 7 8]; % bar x positions
y = [sDurationAvg(plot_order),cDurationAvg(plot_order)]; % bar heights
y_sem = [sDurationSEM(plot_order),cDurationSEM(plot_order)]; % errorbars
max_y = zeros(1,3); % hold highest point in each bar plot

for i = 1:3 % for each state
    face_colors = [1 1 1; colors(plot_order(i),:)];
    for j = 1:2 % for each treatment
        hold on, bar(x(i,j),y(i,j),'EdgeColor',colors(plot_order(i),:),...
            'FaceColor',face_colors(j,:)); % make the barplot
        hold on, errorbar(x(i,j),y(i,j),y_sem(i,j),'LineStyle','none','Color','k') % plot error bars
    end
    % plot datapoints from individual mice
    for j = 1:n_mice
        hold on, plot(x(i,:),[sDuration(plot_order(i),j), cDuration(plot_order(i),j)],...
            'Marker','o','MarkerFaceColor','w','MarkerEdgeColor',[.5 .5 .5],'LineStyle','--',...
            'Color',[.5 .5 .5])
    end
    % get highest point in each bar plot
    max_y(i) = max([cDurationSEM(plot_order(i))+cDurationAvg(plot_order(i)),...
        sDurationSEM(plot_order(i))+sDurationAvg(plot_order(i)),...
        cDuration(plot_order(i),:), sDuration(plot_order(i),:)]);
end
% make some invisible plots for the sake of the legend
hold on, legend_bar{1} = bar(-10,0,'EdgeColor','k','FaceColor','w');
hold on, legend_bar{2} = bar(-10,0,'EdgeColor','k','FaceColor','k');
xlim([-0.2 9.2]) % adjust axis so they're not visible

% calculate how much to offset significance markers
sig_offset = .05*max(max_y);
% plot significance markers
for i = 1:3
    [~,p] = ttest(sDuration(plot_order(i),:),cDuration(plot_order(i),:));
    disp([titles{plot_order(i)},' bout duration p value: ',num2str(p)])
    if p > 0.05 || isnan(p)
        continue
    else
        mark = p_val_to_marker(p);
    end
    hold on, plot(x(i,:), [max_y(i),max_y(i)]+sig_offset, 'k')
    text(mean(x(i,:)), max_y(i)+1.5*sig_offset,mark,'HorizontalAlignment','center','FontSize',14)
end
ylabel('Bout duration (sec)')
set(gca,'XTick',[1.5 4.5 7.5],'XTickLabel',titles(plot_order),'TickDir','out','FontSize',12) % set x labels
box off
% make the legend
legend([legend_bar{1},legend_bar{2}],{'Saline','CNO'},'FontSize',11)

%% plot the detailed bout length bar plots
figure('Color','w','Position',[470 234 1070 268])
x = reshape(1:(nBins*2),2,nBins)'+(0:(nBins-1))'; % get bar plot x values

for i = 1:3
    y = [sHistAvg(plot_order(i),:)',cHistAvg(plot_order(i),:)']; % bar heights
    y_sem = [sHistSEM(plot_order(i),:)',cHistSEM(plot_order(i),:)']; % errorbars
    max_y = zeros(1,nBins); % hold highest point in each bar plot
    bar_handles = cell(1,2);
    subplot(1,3,i); % make a subplot
    
    face_colors = [1 1 1; colors(plot_order(i),:)]; % bar plot face colors
    for j = 1:2 % for each treatment
        hold on, bar_handles{j} = bar(x(:,j),y(:,j),.2,'EdgeColor',colors(plot_order(i),:),...
            'FaceColor',face_colors(j,:)); % make the barplot
        hold on, errorbar(x(:,j),y(:,j),y_sem(:,j),'LineStyle','none','Color','k') % plot error bars
    end
    % plot datapoints from individual mice
    for j = 1:n_mice
        for k = 1:nBins
            hold on, plot(x(k,:),[sHistos(plot_order(i),k,j), cHistos(plot_order(i),k,j)],...
                'Marker','o','MarkerFaceColor','w','MarkerEdgeColor',[.5 .5 .5],'LineStyle','--',...
                'Color',[.5 .5 .5],'MarkerSize',4)
        end
    end
    % fix y limits
    yl = ylim;
    ylim([0 yl(2)]);
    
    % add x labels
    XTickLabels = cell(1,nBins);
    for k = 1:(nBins-1)
        XTickLabels{k} = ['[',num2str(bout_length_ranges(k)),', ',...
            num2str(bout_length_ranges(k+1)),')'];
    end
    XTickLabels{nBins} = ['>= ',num2str(bout_length_ranges(nBins))];
    set(gca,'XTick',mean(x,2),'XTickLabel',XTickLabels,'XTickLabelRotation', 45);
    
    % might put statistical tests here eventually, might not :D
%     % get highest point in each bar plot
%     max_y(i) = max([cCountSEM(plot_order(i))+cCountAvg(plot_order(i)),...
%         sCountSEM(plot_order(i))+sCountAvg(plot_order(i)),...
%         cCount(plot_order(i),:), sCount(plot_order(i),:)]);
    
    % add a legend
    legend([bar_handles{1},bar_handles{2}],{'Saline','CNO'},'Location','northeast','Box','off')

    title(titles{plot_order(i)},'Color',colors(plot_order(i),:))
    xlabel('Bout length (sec)')
    ylabel('Fraction of bouts')
    box off
end

%% bout number
% ok, next we want to find the number of bouts of each state for each mouse
% we will first find it each day, then average across days to get one value
% per mouse
sCount = zeros(3,n_mice); % preallocate
cCount = zeros(3,n_mice);
for i = 1:3
    for j = 1:n_mice
        for k = 1:length(sBouts{i,j})
        	sCount(i,j) = sCount(i,j) + length(sBouts{i,j}{k});
        end
        for k = 1:length(cBouts{i,j})
        	cCount(i,j) = cCount(i,j) + length(cBouts{i,j}{k});
        end
    end
end

% average across mice
sCountAvg = mean(sCount,2);
cCountAvg = mean(cCount,2);

% find variability (SEM)
sCountSEM = std(sCount,[],2)/sqrt(n_mice);
cCountSEM = std(cCount,[],2)/sqrt(n_mice);

%% plot bar plots of count numbers
figure('Color','w','Position',[1307 331 347 268])
x = [1 2; 4 5; 7 8]; % bar x positions
y = [sCountAvg(plot_order),cCountAvg(plot_order)]; % bar heights
y_sem = [sCountSEM(plot_order),cCountSEM(plot_order)]; % errorbars
max_y = zeros(1,3); % hold highest point in each bar plot

for i = 1:3 % for each state
    face_colors = [1 1 1; colors(plot_order(i),:)];
    for j = 1:2 % for each treatment
        hold on, bar(x(i,j),y(i,j),'EdgeColor',colors(plot_order(i),:),...
            'FaceColor',face_colors(j,:)); % make the barplot
        hold on, errorbar(x(i,j),y(i,j),y_sem(i,j),'LineStyle','none','Color','k') % plot error bars
    end
    % plot datapoints from individual mice
    for j = 1:n_mice
        hold on, plot(x(i,:),[sCount(plot_order(i),j), cCount(plot_order(i),j)],...
            'Marker','o','MarkerFaceColor','w','MarkerEdgeColor',[.5 .5 .5],'LineStyle','--',...
            'Color',[.5 .5 .5])
    end
    % get highest point in each bar plot
    max_y(i) = max([cCountSEM(plot_order(i))+cCountAvg(plot_order(i)),...
        sCountSEM(plot_order(i))+sCountAvg(plot_order(i)),...
        cCount(plot_order(i),:), sCount(plot_order(i),:)]);
end
% make some invisible plots for the sake of the legend
hold on, legend_bar{1} = bar(-10,0,'EdgeColor','k','FaceColor','w');
hold on, legend_bar{2} = bar(-10,0,'EdgeColor','k','FaceColor','k');
xlim([-0.2 9.2]) % adjust axis so they're not visible

% calculate how much to offset significance markers
sig_offset = .05*max(max_y);
% plot significance markers
for i = 1:3
    [~,p] = ttest(sCount(plot_order(i),:),cCount(plot_order(i),:));
    disp([titles{plot_order(i)},' bout number p value: ',num2str(p)])
    if p > 0.05 || isnan(p)
        continue
    else
        mark = p_val_to_marker(p);
    end
    hold on, plot(x(i,:), [max_y(i),max_y(i)]+sig_offset, 'k')
    text(mean(x(i,:)), max_y(i)+1.5*sig_offset,mark,'HorizontalAlignment','center','FontSize',14)
end
ylabel(['Number of bouts in ',num2str(bout_end_time),' minutes'])
set(gca,'XTick',[1.5 4.5 7.5],'XTickLabel',titles(plot_order),'TickDir','out','FontSize',12) % set x labels
box off
% make the legend
legend([legend_bar{1},legend_bar{2}],{'Saline','CNO'},'FontSize',11)

%% NREM latency
% convert to minutes
sLatency  = sLatency * epoch_length / 60;
cLatency  = cLatency * epoch_length / 60;

% average across mice
sLatencyAvg = mean(sLatency);
cLatencyAvg = mean(cLatency);

% find variability (SEM)
sLatencySEM = std(sLatency)/sqrt(n_mice);
cLatencySEM = std(cLatency)/sqrt(n_mice);

%% bar plot of NREM latencies
figure('Color','w','Position',[1430 254 255 268])
x = [1 2]; % bar x positions
y = [sLatencyAvg,cLatencyAvg]; % bar heights
y_sem = [sLatencySEM,cLatencySEM]; % errorbars

face_colors = [1 1 1; colors(3,:)];
bar_handles = cell(1,2);
for j = 1:2 % for each treatment
    hold on, bar_handles{j} = bar(x(j),y(j),'EdgeColor',colors(3,:),...
        'FaceColor',face_colors(j,:)); % make the barplot
    hold on, errorbar(x(j),y(j),y_sem(j),'LineStyle','none','Color','k') % plot error bars
end
% plot datapoints from individual mice
for j = 1:n_mice
    hold on, plot(x,[sLatency(j), cLatency(j)],'Marker','o','MarkerFaceColor',...
        'w','MarkerEdgeColor',[.5 .5 .5],'LineStyle','--','Color',[.5 .5 .5])
end
% get highest point in the bar plot
max_y = max([cLatencySEM+cLatencyAvg,sLatencySEM+sLatencyAvg,...
    cLatency, sLatency]);

% calculate how much to offset significance markers
sig_offset = .05*max_y;
% plot significance marker

[~,p] = ttest(sLatency,cLatency);
disp(['NREM latency p value: ',num2str(p)])
if p <= 0.05 && ~isnan(p)
    mark = p_val_to_marker(p);
    hold on, plot(x, [max_y,max_y]+sig_offset, 'k')
    text(mean(x), max_y+1.5*sig_offset,mark,'HorizontalAlignment','center','FontSize',14)
end

ylabel('Latency to first NREM bout (min)')
set(gca,'XTick',[1 2],'XTickLabel',{'Saline','CNO'},'TickDir','out','FontSize',12) % set x labels
box off


function [data, bouts, latency] = loadData(mouselist,nb,tb,h1,h2,min_1st_nrem_epochs)
% Outputs:
% data: average brain state over time. matrix structure: state x time x mouse
%       (we average across days for each mouse)
% bouts: each entry is a cell array, and each entry in those is a cell where
%        each entry contains all bout lengths for one recording.
%        structure: state X mouse. 
% latency: average (over recordings) latency to first NREM period per mouse

% load data
nMice = length(mouselist); % how many mice there are
allMiceFiles = cell(1, nMice); % locations of files
nDays = zeros(1,nMice); % how many days there are for each mouse
% check that the contents of the mouselist are ok, and gather filenames
for m = 1:nMice % for each mouse
    % find sleep stage and laser data
    temp = load(mouselist{m},'dirlist');
    dirlist = temp.dirlist;
    [problemString, fileNames] = AS_checkList(dirlist, {'labels'});
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

% our first output: state X time (binned) X mouse
data = zeros(3,nb/tb,nMice);
% second output: state X mouse. each entry is a cell array, and each entry
% in those is a cell where each entry contains all bout lengths for one
% recording
bouts = cell(3,nMice);
latency = zeros(1,nMice);
% for each mouse
for i = 1:nMice
    % for each day
    % brain state likelihood at each timepoint, each day
    binnedState = zeros(3,nb/tb,nDays(i));
     % preallocate bouts
    for q = 1:3
        bouts{q, i} = cell(1, nDays(i));
    end
    latencies = zeros(1,nDays(i)); % latency in each recording
    for j = 1:nDays(i)
        % load brain state data
        load(allMiceFiles{i}{j}{1},'labels');
        % extract just the brain state for the appropriate length of time
        state = labels(1:nb);
        % get brain state likelihood at each timepoint
        for k = 1:nb/tb
            for p = 1:3
                binnedState(p,k,j) = sum(state(((k-1)*tb+1):(k*tb))==p)/tb;
            end
        end
        % get bouts
        state = labels(h1:h2); % just get brain state for the required time
        if ~isrow(state) % convert to a row
            state = state';
        end
        stateStr = erase(num2str(state),' '); % convert to a string
        for q = 1:3 % for each state
            % get start and end indices of each bout
            [startIndex,endIndex] = regexp(stateStr,[num2str(q),'+']);
            % add to the list of bout lengths
            bouts{q,i}{j} = endIndex - startIndex + 1;
        end
        % get latency to first NREM bout
        startIndex = regexp(stateStr,repmat('3',1,min_1st_nrem_epochs));
        % show a message if none was found
        if isempty(startIndex)
            latencies(j) = length(stateStr);
            disp(['Mouse ',num2str(i),', recording ',num2str(j),...
                ' has no NREM'])
        else
            latencies(j) = startIndex(1) - 1;
        end
    end
    data(:,:,i) = mean(binnedState,3); % average brain states over days
    latency(i) = mean(latencies); % average NREM latency over days
end

function stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
% https://www.mathworks.com/matlabcentral/fileexchange/6874-two-way-repeated-measures-anova
% Aaron Schurger (2005.02.04) Derived from Keppel & Wickens (2004) "Design and Analysis" ch. 18
F1_lvls = unique(F1);
F2_lvls = unique(F2);
Subjs = unique(S);
a = length(F1_lvls); % # of levels in factor 1
b = length(F2_lvls); % # of levels in factor 2
n = length(Subjs); % # of subjects
INDS = cell(a,b,n); % this will hold arrays of indices
CELLS = cell(a,b,n); % this will hold the data for each subject X condition
MEANS = zeros(a,b,n); % this will hold the means for each subj X condition
for i=1:a % F1
    for j=1:b % F2
        for k=1:n % Subjs
            INDS{i,j,k} = find(F1==F1_lvls(i) & F2==F2_lvls(j) & S==Subjs(k));
            CELLS{i,j,k} = Y(INDS{i,j,k});
            MEANS(i,j,k) = mean(CELLS{i,j,k});
        end
    end
end
AB = reshape(sum(MEANS,3),a,b); % across subjects
AS = reshape(sum(MEANS,2),a,n); % across factor 2
BS = reshape(sum(MEANS,1),b,n); % across factor 1
A = sum(AB,2); % sum across columns, so result is ax1 column vector
B = sum(AB,1); % sum across rows, so result is 1xb row vector
S = sum(AS,1); % sum across columns, so result is 1xs row vector
T = sum(sum(A)); % could sum either A or B or S, choice is arbitrary
dfA = a-1;
dfB = b-1;
dfAB = (a-1)*(b-1);
% dfS = n-1;
dfAS = (a-1)*(n-1);
dfBS = (b-1)*(n-1);
dfABS = (a-1)*(b-1)*(n-1);
expA = sum(A.^2)./(b*n);
expB = sum(B.^2)./(a*n);
expAB = sum(sum(AB.^2))./n;
expS = sum(S.^2)./(a*b);
expAS = sum(sum(AS.^2))./b;
expBS = sum(sum(BS.^2))./a;
expY = sum(sum(sum(MEANS.^2))); %sum(Y.^2);
expT = T^2 / (a*b*n);
ssA = expA - expT;
ssB = expB - expT;
ssAB = expAB - expA - expB + expT;
% ssS = expS - expT;
ssAS = expAS - expA - expS + expT;
ssBS = expBS - expB - expS + expT;
ssABS = expY - expAB - expAS - expBS + expA + expB + expS - expT;
% ssTot = expY - expT;
msA = ssA / dfA;
msB = ssB / dfB;
msAB = ssAB / dfAB;
% msS = ssS / dfS;
msAS = ssAS / dfAS;
msBS = ssBS / dfBS;
msABS = ssABS / dfABS;
fA = msA / msAS;
fB = msB / msBS;
fAB = msAB / msABS;
pA = 1-fcdf(fA,dfA,dfAS);
pB = 1-fcdf(fB,dfB,dfBS);
pAB = 1-fcdf(fAB,dfAB,dfABS);
stats = {'Source','SS','df','MS','F','p';...
    FACTNAMES{1}, ssA, dfA, msA, fA, pA;...
    FACTNAMES{2}, ssB, dfB, msB, fB, pB;...
    [FACTNAMES{1} ' x ' FACTNAMES{2}], ssAB, dfAB, msAB, fAB, pAB;...
    [FACTNAMES{1} ' x Subj'], ssAS, dfAS, msAS, [], [];...
    [FACTNAMES{2} ' x Subj'], ssBS, dfBS, msBS, [], [];...
    [FACTNAMES{1} ' x ' FACTNAMES{2} ' x Subj'], ssABS, dfABS, msABS, [], []};

function [mark] = p_val_to_marker(p)
% decide which marker to plot for different p values
if p <= .05
    mark = '*';
end
if p < .01
    mark = '**';
end
if p < .001
    mark = '***';
end
