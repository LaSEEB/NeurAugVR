clear, clc, close all
% The purpose of this script is to show how the ERSP calculated by
% normalizing each TF-trial and then averaging (instead of the reverse)
% overstimates the synchronization (produces an artificial ERS%)
%% Create EEG
EEG = eeg_emptyset;

fs = 250;
base_dur = 5;
task_dur = 5;
trials_num = 1000;

% Data
mu = 0;
sigma = 100;
signal = normrnd(0,100,1,fs*(base_dur+task_dur)*trials_num);
EEG.data = signal;
EEG.srate = fs;
EEG.times = (0:numel(signal)-1)/fs;

% Events
t = 0;
c = 1;
for i = 1:trials_num
    EEG.event(c).type = 'base';
    EEG.event(c).latency = t * fs;
    EEG.event(c).urevent = c;
    t = t + base_dur;
    c = c + 1;
    EEG.event(c).type = 'task';
    EEG.event(c).latency = t * fs;
    EEG.event(c).urevent = c;
    t = t + task_dur;
    c = c + 1;
end

%% Calculate + plot ERSP
% Epoch
dire = {'task'};  % Left
epochlims = [-5, 5];
EEGdi = pop_epoch(EEG, dire, epochlims, 'epochinfo', 'yes');
% Select channel
% chan = 'C4';
% chn = find(ismember({EEG.chanlocs(:).labels}, chan));
chn = 1;
tubmat = double(EEGdi.data(chn, :, :));  % I use the word 'tub' for a 1x1xN array, and 'tubmat' for a 1xMxN array
% Time-Frequency decomposition of each trial
baseline = nan; % Indiferent actually
scale = 'abs';
cycles = 0;
[~, ~, ~, time_vec, freq_vec, ~, ~, amp_vol] = newtimef_trueamp(tubmat, size(tubmat,2), [EEGdi.times(1) EEGdi.times(end)], EEGdi.srate, cycles,'plotitc','off', 'baseline', baseline, 'scale', scale, 'plotersp', 'off');
time_vec = time_vec/1000;
amp_vol = double(amp_vol);
% Calculate power
pow_vol = amp_vol.^2;
% Average trials
pow_mat = mean(pow_vol,3);
% Calculate ERD%
baselims = [-4, 0];
erd_mat = (pow_mat./mean(pow_mat(:,time_vec > baselims(1) & time_vec < baselims(2)),2) - 1)*100;
% Average alpha's frequency band
freqlims = [8, 12];
erd_vec = mean(erd_mat(freq_vec > freqlims(1) & freq_vec < freqlims(2),:),1);
% Average task's time interval
tasklims = [0.5, 4];
erd_val = mean(erd_vec(time_vec > tasklims(1) & time_vec < tasklims(2)));
% Plot ERSP%
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

%% Calculate + plot ERSP [normalization: single-trial]
erd_vol = (pow_vol./mean(pow_vol(:,time_vec > baselims(1) & time_vec < baselims(2),:),2) - 1)*100;
erd_mat_s = mean(erd_vol,3);
figure
imagesc('XData',time_vec,'YData',freq_vec,'CData',erd_mat_s,clims)
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

%%
% pow_bases = mean(mean(pow_vol(:,time_vec > baselims(1) & time_vec < baselims(2),:),2),1);
% pow_tasks = mean(mean(pow_vol(:,time_vec > tasklims(1) & time_vec < tasklims(2),:),2),1);
% figure
% histogram(pow_bases(:))
% hold on
% histogram(pow_tasks(:))
%%
% erd_vals = (pow_tasks(:)-pow_bases(:))./pow_bases(:)*100;
% erd_vals = mean(mean(erd_vol(:,time_vec > tasklims(1) & time_vec < tasklims(2),:),2),1);
% figure
% histogram(erd_vals(:))