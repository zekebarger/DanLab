% zeke barger 010720
% Compare sleeping behavior between control and treatment groups over
% long time scales (for example, chemogenetic manipulations)
%
% NOTE: the statistics are not quite ready for prime time, since
% there's no correction for multiple comparisons--just a
% paired t-test at each timepoint. 
% For the bout length histogram, bouts are pooled across days to make a
% histogram for each mouse. Then the histograms are averaged across mice.
% Right now, the 'histogram' is shown as a line plot. This can be easily
% changed - overlapping transparent histograms might be preferable.
%

function [] = AS_chemoPlot()
%% user-defined parameters
hrs2take = 3; % hours of data to take from the beginning of each recording
epochLength = 2.5; % length of brain state epochs, in seconds
timebin = 20; % minutes of brain states to avg for each plotted timepoint
controlLabel = 'Saline'; % label for control experiments (e.g., 'saline')
treatmentLabel = 'CNO'; % label for treatment experiments (e.g., 'CNO')
plotOrder = [3 2 1]; % order to plot the brain states
                     % 1 = REM, 2 = wake, 3 = NREM
% the two parameters below are the times, in minutes, between which to
% collect data for the histogram of bout durations (a bout is an 
% uninterrupted period of some brain state)
histoStart = 0; % time to start, in minutes
histoEnd = 180; % time to end, in minutes
hBinWidth = 10; % bin width for histogram, in seconds
hMax = 3*60; % maximum bout length in the histogram, in seconds

%% extract data from recordings
% convert from hours/minutes to epochs
nb = round(hrs2take*60*60/epochLength); % epochs from beginning
tb = timebin*60/epochLength; % epochs per timepoint
h1 = max([round(histoStart*60/epochLength), 1]); % epoch to start histogram
h2 = round(histoEnd*60/epochLength); % epoch to stop histogram


% have the user select the mouselists for each type of experiment
[sFileName,sPathName,~] = uigetfile('*.mat','Choose CONTROL (e.g., saline) mouselist');
[cFileName,cPathName,~] = uigetfile([sPathName, '*.mat'],...
    'Choose TREATMENT (e.g., CNO) mouselist');

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

% load brain state data and bout lengths
[sData,sBouts] = loadData(sList.mouselist,nb,tb,h1,h2);
[cData,cBouts] = loadData(cList.mouselist,nb,tb,h1,h2);
sData = sData*100;
cData = cData*100;

% if there was some problem
if isempty(sData) || isempty(cData)
    return
end

% average across mice
sAvg = mean(sData,3);
cAvg = mean(cData,3);

% find variability (SEM)
sSEM = std(sData,[],3)/sqrt(size(sData,3));
cSEM = std(cData,[],3)/sqrt(size(cData,3));

% get bout length histograms
edges = 0:hBinWidth:hMax; % compute histogram edges
nBins = length(edges) - 1; % get number of bins
centers = edges(1:end-1) + hBinWidth/2; % find centers of bins
% preallocate for the histograms from each mouse
sHistos = zeros(3,nBins,size(sBouts,2));
cHistos = zeros(3,nBins,size(cBouts,2));
% calculate histogram for each mouse
for i = 1:3 % for each state
    for j = 1:size(sBouts,2)
        [sHistos(i,:,j),~] = histcounts(sBouts{i,j}, edges,'Normalization','probability');
    end
    for j = 1:size(cBouts,2)
        [cHistos(i,:,j),~] = histcounts(cBouts{i,j}, edges,'Normalization','probability');
    end
end

% average across mice
sHistAvg = mean(sHistos,3);
cHistAvg = mean(cHistos,3);

% find variability (SEM)
sHistSEM = std(sHistos,[],3)/sqrt(size(sBouts,3));
cHistSEM = std(cHistos,[],3)/sqrt(size(cBouts,3));

%% plot the brain state likelihoods
% set up time axis
t = (timebin * (1:(nb/tb)) +5 - timebin/2) / 60 ;

titles = {'REM','Wake','NREM'};
colors = [43 160 220;
    129 131 132
    255 194 61] ./ 255; % rem wake nrem

figure('Color','w','Position',[458 408 1010 240])
plotCount = 1; % keeps track of how many plots we've made so far
plotHandles = [0 0 0]; % hold handles to our plots
for i = plotOrder
    plotHandles(i) = subplot(1,3,plotCount); % make a subplot
    
    % plot mean and sem
    errorbar(t,sAvg(i,:),sSEM(i,:),'o-','Color',colors(i,:),...
        'CapSize',0,'MarkerFaceColor','w',...
        'LineWidth',1);
    hold on, errorbar(t,cAvg(i,:),cSEM(i,:),'o-','Color',...
        colors(i,:),'MarkerFaceColor',colors(i,:),...
        'LineWidth',1,'CapSize',0);
    
    % show significance for each timepoint
    for x = 1:length(t)
        [~,p] = ttest(squeeze(sData(i,x,:)),squeeze(cData(i,x,:)));
        
        if p > 0.05 || isnan(p)
            continue
        else
            if p<=.05
                mark = '*';
            end
            if p<.01
                mark = '**';
            end
            if p < .001
                mark = '***';
            end
        end
        yl = ylim;
        y=max([sAvg(i,x)+sSEM(i,x), cAvg(i,x)+cSEM(i,x)])+.03*diff(yl);
        
        text(t(x), y,mark,'HorizontalAlignment','center','FontSize',12)
        
    end
    
    % add a legend
    if plotCount==1
        legend({controlLabel,treatmentLabel},'Location','northeast','Box','off')
        ylabel('Percentage (%)')
    end
    title(titles{i},'Color',colors(i,:))
    xlabel('Time (h)')
    box off
    
    % set lower y limit to 0
    yl = ylim;
    set(gca,'YLim',[0 yl(2)]);
    
    % set y axis for wake and nrem to be the same
    if plotCount == 3 % if both wake and nrem have been plotted
        % get existing y limits
        wY = ylim(plotHandles(2));
        nY = ylim(plotHandles(3));
        % get max value for upper limits
        yNew = [0 max([wY(2) nY(2)])];
        % apply new values
        ylim(plotHandles(2), yNew);
        ylim(plotHandles(3), yNew);
    end
    plotCount = plotCount + 1; % update number of plots
end

%% plot the bout length histograms
figure('Color','w','Position',[458 90 1010 240])
plotCount = 1; % keeps track of how many plots we've made so far
plotHandles = [0 0 0]; % hold handles to our plots
for i = plotOrder
    plotHandles(i) = subplot(1,3,plotCount); % make a subplot
    
    % plot mean and sem
    errorbar(centers,sHistAvg(i,:),sHistSEM(i,:),'o-','Color',colors(i,:),...
        'CapSize',0,'MarkerFaceColor','w',...
        'LineWidth',1);
    hold on, errorbar(centers,cHistAvg(i,:),cHistSEM(i,:),'o-','Color',...
        colors(i,:),'MarkerFaceColor',colors(i,:),...
        'LineWidth',1,'CapSize',0);
    
    % show significance for each timepoint
    for x = 1:length(centers)
        [~,p] = ttest(squeeze(sHistos(i,x,:)),squeeze(cHistos(i,x,:)));
        
        if p > 0.05 || isnan(p)
            continue
        else
            if p<=.05
                mark = '*';
            end
            if p<.01
                mark = '**';
            end
            if p < .001
                mark = '***';
            end
        end
        yl = ylim;
        y=max([sHistAvg(i,x)+sHistSEM(i,x), cHistAvg(i,x)+cHistSEM(i,x)])+.03*diff(yl);
        
        text(centers(x), y,mark,'HorizontalAlignment','center','FontSize',12)
        
    end
    
    % add a legend
    if plotCount==1
        legend({controlLabel,treatmentLabel},'Location','northeast','Box','off')
        ylabel('Percentage (%)')
    end
    title(titles{i},'Color',colors(i,:))
    xlabel('Bout length (sec)')
    ylabel('Fraction of bouts')
    box off
    
    plotCount = plotCount + 1; % update number of plots
end


%% function to load data from files
% return a matrix: state x time x mouse
% average across days for each mouse
function [brainStateData, bouts] = loadData(mouselist,nb,tb,h1,h2)
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
brainStateData = zeros(3,nb/tb,nMice);
% second output: state X mouse
bouts = cell(3,nMice);
% for each mouse
for i = 1:nMice
    % for each day
    % brain state likelihood, state X time (binned) X recording
    binnedState = zeros(3,nb/tb,nDays(i));
    % preallocate bouts
    for q = 1:3
        bouts{q, i} = [];
    end
    for j = 1:nDays(i)
        % load brain state data
        load(allMiceFiles{i}{j}{1},'labels');
        % check if the recording is long enough
        if length(labels) < nb || length(labels) < h2
           disp(['Error: ',allMiceFiles{i}{j}{1},' is too short.']) 
           brainStateData = [];
           return
        end
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
            bouts{q,i} = [bouts{q,i}, endIndex - startIndex + 1];
        end
    end
    % average across recordings
    brainStateData(:,:,i) = mean(binnedState,3);
end
