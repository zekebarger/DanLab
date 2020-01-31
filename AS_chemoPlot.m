% Zeke Barger, 013120
% plot brain states following chemogenetic manipulation
% NOTE: the saline and CNO mouselists should point to the same mice,
% in the same order!
% Significance is calculated using a two-way repeated measures ANOVA
% with Bonferroni correction followed by paired t-tests
% *: p < 0.05, **: p < 0.01, ***: p < 0.001
% It's assumed that each recording starts right after the injection

function [] = AS_chemoPlot()
%% user-defined parameters
% how many hrs of data to take from the beginning of each recording
hours_to_plot = 3;
epoch_length = 2.5; % in seconds
mins_per_bin = 20; % minutes of brain states to avg for each plotted timepoint
plot_order = [3 2 1]; % 1=R, 2=W, 3=N. So, to plot wake, then NREM, then REM,
%                       set this to [2 3 1]
recording_delay = 5; % minutes between injection and recording start

%% get brain state data
% number of bins to take - don't change this
total_epochs = round(hours_to_plot*60*60/epoch_length); % bins from beginning
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

% load brain state data
sData = loadData(sList.mouselist,total_epochs,epochs_per_bin)*100;
cData = loadData(cList.mouselist,total_epochs,epochs_per_bin)*100;

% if there was some problem
if isempty(sData) || isempty(cData)
    return
end

% average across mice
sAvg = mean(sData,3);
cAvg = mean(cData,3);

% find variability
sSEM = std(sData,[],3)/sqrt(size(sData,3));
cSEM = std(cData,[],3)/sqrt(size(cData,3));

%% plot
% set up time axis
t = (mins_per_bin * (1:(total_epochs/epochs_per_bin)) + recording_delay - mins_per_bin/2) / 60 ;

titles = {'REM','Wake','NREM'};
colors = [43 160 220;
    129 131 132
    255 194 61] ./ 255; % rem wake nrem
%%
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
        all_obj = {};
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
    disp([titles{i},': F(',df1,',',df2,') = ',F,', p = ',pval_string])
    % if there was a significant treatment x time interaction, run post-hoc
    % tests at each timepoint
    if pval <= (0.05 / 3) % Bonferroni correction
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
            y=max([sAvg(i,x)+sSEM(i,x), cAvg(i,x)+cSEM(i,x)])+.028*diff(yl);
            text(t(x), y,mark,'HorizontalAlignment','center','FontSize',14)         
        end
    end 
    plot_number = plot_number + 1;
end
%%
% return a matrix: state x time x mouse
% average across days for each mouse
function [data] = loadData(mouselist,nb,tb)
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

data = zeros(3,nb/tb,length(mouselist));
% for each mouse
for i = 1:nMice
    % for each day
    % brain state likelihood at each timepoint, each day
    binnedState = zeros(3,nb/tb,nDays(i));
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
    end
    
    data(:,:,i) = mean(binnedState,3);
end

function stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
% https://www.mathworks.com/matlabcentral/fileexchange/6874-two-way-repeated-measures-anova
% Aaron Schurger (2005.02.04) Derived from Keppel & Wickens (2004) "Design and Analysis" ch. 18
stats = cell(4,5);
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
dfS = n-1;
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
ssS = expS - expT;
ssAS = expAS - expA - expS + expT;
ssBS = expBS - expB - expS + expT;
ssABS = expY - expAB - expAS - expBS + expA + expB + expS - expT;
ssTot = expY - expT;
msA = ssA / dfA;
msB = ssB / dfB;
msAB = ssAB / dfAB;
msS = ssS / dfS;
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
