% zeke barger 121319
% plots some example data during laser stimulation

function AS_exampleLaserPlot()
%% parameters - you should modify these
dpath = 'D:\Data\other\AStestdata\Dex1\2018-07-30_11-27-09'; % folder containing the data
t1 = 5; % time to start plotting, in minutes
t2 = 120; % time to stop plotting, in minutes
epochLength = 2.5; % brain state epoch length
SR = 512; % EEG/EMG sampling rate
f1 = 1; % first frequency band of spectrogram to plot, in Hz
f2 = 25;
mVScaleBar = 2; % length of EMG scale bar in mV
mV = 1000; % for EMG, the value equivalent to 1 mV
% brain state colors (REM, wake, NREM)
colors = [43 160 220;
            129 131 132
            255 194 61
            0 0 0] ./ 255;

%% process data
% get names of files
[problemString, fileNames] = AS_checkEntry(dpath, {'EEG','EMG','laser','labels'});

% if something is missing, quit
if ~isempty(problemString)
    disp(problemString)
    return
end

% load the data
temp = load(fileNames{1},'EEG');
EEG = temp.EEG;
temp = load(fileNames{2},'EMG');
EMG = temp.EMG;
temp = load(fileNames{3},'laser');
laser = temp.laser;
temp = load(fileNames{4},'labels');
labels = temp.labels;

% calculate spectrogram
[s, t, f] = createSpectrogram(EEG, SR, epochLength);

% get index of t1, t2 epochs
t1e = floor(t1*60/epochLength)+1;
t2e = floor(t2*60/epochLength)+1;

% get index of t1, t2 in the EEG/EMG data
t1m = round(((t1e-1)*epochLength)*SR);
t2m = round(t2e*epochLength*SR);

% trim EMG
EMG = EMG(t1m:t2m);

% trim spectrogram and labels
[~,f1Idx] = min(abs(f - f1)); 
[~,f2Idx] = min(abs(f - f2)); 
s = s(t1e:t2e, f1Idx:f2Idx);
labels = labels(t1e:t2e);

% make time axis, in epochs
t = round(t*10)/10;
t = t(t1e:t2e);
te = t - t(1) + epochLength/2;

% make time axis, in seconds
ts = ((1:length(EMG)) - 1) / SR;

% get laser onsets and offsets
i = 1;
while i <= size(laser,1) % for each stimulus
    % if offset is before t1, or onset after t2, discard
    if laser(i,2) < t1*60 || laser(i,1) > t2*60 
        laser(i,:) = [];
        continue
    end
    % if onset is before t1, set to t1
    if laser(i,1) < t1*60 
        laser(i,1) = t1;
    end
    % if offset is after t2, set to t2
    if laser(i,2) > t2*60 
        laser(i,2) = t2;
    end 
    i = i+1;
end
if ~isempty(laser)
    laser = laser - (t1e-1)*epochLength;
end

%% plot
figure('Color','w','Position',[150,370,1600,420])

% plot laser
if ~isempty(laser)
    axes('Position',[.05 .05 .9 .9]);
    for i = 1:size(laser,1)
        hold on, fill([laser(i,1),laser(i,1),laser(i,2),laser(i,2)],[0 .96 .96 0],...
            [.6 .59 1],'EdgeAlpha',0);
        if i == 1
            text(mean(laser(1,:)), .99,'Laser','HorizontalAlignment','center',...
                'FontSize',14,'Color',[.6 .59 1])
        end
    end
    xlim([ts(1), ts(end)])
    ylim([0 1])
    set(gca,'XTick',[],'YTick',[]);
    axis off
end

% spec axis
axes('Position',[.05 .42 .9 .45]);
set(gca,'YTick',[f1 f2],'XTick',[],'FontSize',14,'TickLength',[0 0])
ylim([f1 f2])
ylabel('Freq (Hz)')

% plot spectrogram
pLow = prctile(s,1,'all');
High = prctile(s,99,'all');
axes('Position',[.05 .42 .9 .45]);
imagesc(te,f,s'); axis xy
set(gca,'XTick',[],'YTick',[],'CLim',[pLow,High]); 
axis off
colorbar('Location','north','Position',[0.9,0.88,0.05,0.03],'Ticks',[],'Color','none')
text(.005*te(end),.92*f(end),'EEG','Color','w','FontSize',14)

% plot emg
axes('Position',[.05 .16 .9 .24]);
plot(ts,EMG,'k')
xlim([ts(1), ts(end)])
set(gca,'XTick',[],'YTick',[]); 
axis off
yl = ylim;
xl = xlim;
text(-.04*ts(end), 0, 'EMG','FontSize',14,'VerticalAlignment','middle')

% EMG scale bar
hold on, plot(ones(1,2)*.001*(xl(2)-xl(1)), [yl(1), yl(1)+mVScaleBar*mV],'k','LineWidth',2)
text(.007*(xl(2)-xl(1)), yl(1)+mVScaleBar*mV/2, [num2str(mVScaleBar), ' mV'],...
    'FontSize',14)

% time scale bar
hold on, plot([0 60], ones(1,2)*.98*yl(2),'k','LineWidth',2)
text(0, .7*yl(2), '60 s','FontSize',14)

% plot brain state
aa = axes('Position',[.05 .07 .9 .07]);
for i = 3:-1:1
    hold on, plot(0,0,'Color',colors(i,:),'LineWidth',5)
end
hold on, imagesc(te,1,labels,[1 4]);
colormap(aa,colors);
set(gca,'XTick',[],'YTick',[]); 
axis off
xlim([ts(1), ts(end)])
ylim([.9, 1.1])
legend({'NREM','Wake','REM'},'Orientation','horizontal','Position',[0.53,0.013,0.41,0.04],...
    'FontSize',14)
legend('boxoff')
