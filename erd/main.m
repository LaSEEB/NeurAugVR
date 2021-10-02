%% Add paths
clear, clc
addpath('C:/Users/guta_/Documents/MATLAB/eeglab2021.0')
varsbefore = who; eeglab; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})

%% Load
EEG = pop_loadbv('..\data\', 'sub-pilot03_ses-outside_task-neurowMIMO.vhdr', [], []);

%% Preprocess (prep. 1)
bandlims = [1, 40];
resamp = 250;
EEG = pop_select(EEG, 'nochannel',{'ECG','EKG'});
EEG = pop_resample(EEG, resamp);
EEG = pop_eegfiltnew(EEG, 'locutoff',bandlims(1), 'plotfreqz',0);
EEG = pop_eegfiltnew(EEG, 'hicutoff',bandlims(2), 'plotfreqz',0);
EEG = pop_reref(EEG, []);

%% Epoch
dire = {'S  7'};  % Left
epochlims = [-5.5, 5.5];
EEGdi = pop_epoch(EEG, dire, epochlims, 'epochinfo', 'yes');

%% Select channel
chan = 'C4';
chn = find(ismember({EEG.chanlocs(:).labels}, chan));
tubmat = double(EEGdi.data(chn, :, :));  % I use the word 'tub' for a 1x1xN array, and 'tubmat' for a 1xMxN array

%% Time-Frequency decomposition of each trial
baseline = nan; % Indiferent actually
scale = 'abs';
cycles = 0;
[~, ~, ~, time_vec, freq_vec, ~, ~, amp_vol] = newtimef_trueamp(tubmat, size(tubmat,2), [EEGdi.times(1) EEGdi.times(end)], EEGdi.srate, cycles,'plotitc','off', 'baseline', baseline, 'scale', scale, 'plotersp', 'off');
time_vec = time_vec/1000;
amp_vol = double(amp_vol);

%% Calculate power
pow_vol = amp_vol.^2;

%% Average trials
pow_mat = mean(pow_vol,3);

%% Calculate ERD%
baselims = [-4, 0];
erd_mat = (pow_mat./mean(pow_mat(:,time_vec > baselims(1) & time_vec < baselims(2)),2) - 1)*100;

%% Average alpha's frequency band
freqlims = [8, 12];
erd_vec = mean(erd_mat(freq_vec > freqlims(1) & freq_vec < freqlims(2),:),1);

%% Average task's time interval
tasklims = [0.5, 4];
erd_val = mean(erd_vec(time_vec > tasklims(1) & time_vec < tasklims(2)));

%% Plot ERSP%
xlims = [time_vec(1), time_vec(end)];
ylims = [0, 30];
clims = [-100,100];
figure
imagesc('XData',time_vec,'YData',freq_vec,'CData',erd_mat,clims)
hold on
colormap jet
ch = colorbar;
set(get(ch,'title'),'string','ERSP%');
% caxis(clims)
plot(xlims,freqlims(1)*[1,1],'Color',[0,0,0])
hold on
plot(xlims,freqlims(2)*[1,1],'Color',[0,0,0])
hold on
plot([0,0],ylims,'Color',[0,0,0],'LineStyle','--')
hold on
ylim(ylims)
xlim(xlims)
xlabel('Time [s]')
ylabel('Frequency [Hz]')

%% Plot time-Power
pow_vec = mean(pow_mat(freq_vec > freqlims(1) & freq_vec < freqlims(2),:),1);
pow_std_vec = std(pow_mat(freq_vec > freqlims(1) & freq_vec < freqlims(2),:),0,1);
pow_val = mean(pow_vec(time_vec > tasklims(1) & time_vec < tasklims(2)));
top_margin_vec = pow_vec + pow_std_vec;
bot_margin_vec = pow_vec - pow_std_vec;

colors = lines(2);
figure
p(1) = plot(time_vec, pow_vec, 'Color', colors(1,:));  % use colors(2,:) for orange color
hold on
h = fill([time_vec, time_vec(end:-1:1)], [bot_margin_vec, top_margin_vec(end:-1:1)], colors(1,:),'LineStyle','none'); % we make the shape to fill starting w/ bottom edge from left to right, and then top edge from right to left
set(h,'facealpha',.3)
hold on
text(time_vec(end),pow_val,sprintf(' %d µV^2', round(pow_val)),'FontSize',10, 'Color', colors(1,:))
xlim(xlims)
xlabel('Time [s]')
ylabel('Power [µV^2]')
legend([p(1)],'Left - C4')




