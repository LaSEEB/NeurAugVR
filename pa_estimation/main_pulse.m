clc, clear, close all
%% Add paths
rfol = 'C:/Users/guta_/Desktop/Data Analysis/';             % Root folder
cfol = fileparts(matlab.desktop.editor.getActiveFilename);  % Current file folder
restoredefaultpath; clear RESTOREDEFAULTPATH_EXECUTED;      % Remove any added paths
cd(cfol);                                                   % Change directory to the executing file's directory
addpath(strcat(rfol,'Libraries/eeglab2021.0'))              % Add EEGLAB            
varsbefore = who; eeglab; close; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:});  % Run EEGLAB but remove variables and close figure
clear rfol cfol

%% Load data
load('sub-16_ses-inside2_task-neurowMIMO_run-01_qrs','EEG');
load('sub-16_ses-inside2_task-neurowMIMO_run-01_puls','puls');
EEG = pop_eegfiltnew(EEG, 'locutoff',0.5, 'plotfreqz',0);
EEG = pop_eegfiltnew(EEG, 'hicutoff',40, 'plotfreqz',0);

%% Extract
% Signals
eeg = EEG.data(ismember(upper({EEG.chanlocs(:).labels}),{'T7'}),:);
ecg = EEG.data(ismember(upper({EEG.chanlocs(:).labels}),{'ECG','EKG'}),:);
ppg = puls.data;
fs = EEG.srate;

% Markers
ecg_markers = double([EEG.event(strcmp('QRSi', {EEG.event(:).type})).latency]);
ppg_markers = puls.peaksmarkers;

%% Match EEG, ECG and PPG
% The EEG and ECG are sampled with the same sampling rate and alligned already
% The PPG is neither. Also, its markers don't hit its maxima exactly. So:
% 1: resample;
% 2: snap them to the maxima;
% 3: allign!

% Resample
eeg_time = (0:numel(eeg)-1)/fs;
ppg1 = ppg;
ppg1_time = (0:numel(ppg)-1)/puls.fs;
ppg1_marker = ppg1(puls.peaksmarkers);
ppg1_marker_time = ppg1_time(puls.peaksmarkers);
ppg2 = resample(ppg,fs,puls.fs);
ppg2_time = (0:numel(ppg2)-1)/fs;
ppg2_marker_ids = dsearchn(ppg2_time',ppg1_marker_time');
ppg2_marker = ppg2(ppg2_marker_ids);
ppg2_marker_time = ppg2_time(ppg2_marker_ids);


% Snap markers to maxima
L_data = numel(ppg2_time);
snap_nhood = 0.1; % seconds, ~2%
snap_margins = ceil(snap_nhood * fs / 2) * ones(1,numel(ppg2_marker_ids));
low_mask = (ppg2_marker_ids-snap_margins)<1;
high_mask = (ppg2_marker_ids+snap_margins)>L_data;
snap_margins(low_mask) = ppg2_marker_ids(low_mask)-1;
snap_margins(high_mask) = L_data - ppg2_marker_ids(high_mask);
ppg3_marker_ids = zeros(size(ppg2_marker_ids));
for i = 1:numel(ppg2_marker_ids)
    lat = ppg2_marker_ids(i);
    snap_margin = snap_margins(i);
    [~,I] = max(ppg2(lat - snap_margin : lat + snap_margin));
    ppg3_marker_ids(i) = lat - snap_margin - 1 + I;
end

ppg3_marker = ppg2(ppg3_marker_ids);
ppg3_marker_time = ppg2_time(ppg3_marker_ids);

% Plot the PPG transformations so far
figure
subplot(3,1,1)
plot(ppg1_time,ppg1)
hold on
plot(ppg1_marker_time,ppg1_marker,'r*')
title('Original PPG + markers')
subplot(3,1,2)
plot(ppg2_time,ppg2)
hold on
plot(ppg2_marker_time,ppg2_marker,'r*')
title('Resampled PPG + markers')
subplot(3,1,3)
plot(ppg2_time,ppg2)
hold on
plot(ppg3_marker_time,ppg3_marker,'r*')
title('Resampled PPG + snapped markers')
xlabel('Time [s]')

% Allign (and trim) PPG to EEG(/ECG)
ppg4_time = ppg2_time - puls.delay;
ppg_eeg_mask = ppg4_time>=0 & ppg4_time<=eeg_time(end);
ppg4_time = ppg4_time(ppg_eeg_mask);
ppg4 = ppg2(ppg_eeg_mask);

ppg3_marker_ids_logical = zeros(1,numel(ppg2_time));
ppg3_marker_ids_logical(ppg3_marker_ids) = 1;
ppg4_markers_logical = ppg3_marker_ids_logical(ppg_eeg_mask);
ppg4_marker_ids = find(ppg4_markers_logical);
ppg4_marker = ppg4(ppg4_marker_ids);
ppg4_marker_time = ppg4_time(ppg4_marker_ids);

% Plot EEG, ECG and PPG alligned
figure
subplot(3,1,1)
plot(eeg_time,eeg)
title('EEG')
subplot(3,1,2)
plot(eeg_time,ecg)
title('ECG')
hold on
plot(eeg_time(ecg_markers),ecg(ecg_markers),'r*')
subplot(3,1,3)
plot(ppg4_time,ppg4)
hold on
plot(ppg4_marker_time,ppg4_marker,'r*')
title('PPG')
xlabel('Time [s]')

%% Evaluate how well the ECG estimates the PA on the EEG
% The idea is to use the ECG R peaks to segment the EEG and compute a
% metric which tells how well those segments estimate the PA. 
% (then to do the same with the PPG peaks, and compare the metrics, to conclude if it's better to use the ECG or the PPG)
% Below is the computation of the metric using the ECG R peaks. It can and
% should be tweaked, specially the choice of metric. The computation using
% the PPG is not done, but the code should be the same, just using the PPG
% markers instead of the ECG R peaks
% So to not pre-determine where to segment the PA segments given the R
% peaks (from R to R? from R + 1 to R + 1? etc.), a set of time shifts
% (from the R peaks) is tried.
% Also, (for each time shift), a number of segments (win_num) is
% extracted repeatedly until the end of the recorder. The metric is
% calculated for each time shift, for each extraction. In the end, to
% compare timeshifts, the median of the metrics are used (boxplots)
%
% Future work: Find a better suited metric and use it to compare how well the ECG
% and PPG catch the PA. Also, a lot of code should be put into single-line
% functions to make it easier to read (like snap_markers_to_maxima(), ...).
% Also, one should make sure the ECG and PPG have all the markers. I think
% this particular PPG is missing 2.
%% Set parameters
time_shift_durs = -0.5:0.1:0.5;  % Time shifts from the R peak to try when segmenting the EEG (e.g. 0 => segment goes from one R peak to the next)
win_max_dur = 1.5;  % Max cardiac cycle duration (and hence max template duration)
win_num = 20;  % Number of segments for one template

%% Segment PA
% For each time shift (in time_shift_durs), 
% for each position (in slides_num) that the segments (win_num) can be
% extracted: extract the segments, average them into a template and compute a
% metric to evaluate how well they catch the PA
time_shift_lens = round(time_shift_durs*fs);
slides_num = numel(ecg_markers)-win_num;
win_max_len = win_max_dur*fs;
metric = zeros(numel(time_shift_lens), slides_num);
avg_wins = zeros(numel(time_shift_lens), win_max_len, slides_num);

for i = 1:numel(time_shift_lens)
    for j = 1:slides_num
        ecg_markers_win = ecg_markers(j:j+win_num) + time_shift_lens(i);
        ecg_markers_win = ecg_markers_win(ecg_markers_win>0 & ecg_markers_win<numel(eeg));
        wins = segment_signal(eeg,ecg_markers_win,win_max_len);
        avg_win = mean(wins,1);
        metric(i,j) = mean((wins).^2,'all','omitnan');  % Compute a quality metric for each. In this case it's the average squared segments
        avg_wins(i,:,j) = avg_win;
    end
end

%% Plot PA segments
times = (0:(win_max_len-1))/fs;
figure
for i = 1:numel(time_shift_lens)
    subplot(1,numel(time_shift_lens),i)
    for j = 1:slides_num
           plot(times,avg_wins(i,:,j))
           hold on
    end
end

%% Plot metrics
% For each time shift (in time_shift_durs), 
% for each position (in slides_num) at which the segments (win_num) are extracted, 
% plot a metric evaluating how well the segments catch the PA
figure
boxmat = metric';
boxplot(boxmat,'colors',[0.5,0.5,0.5],'labels',time_shift_durs)
hold on
plot(median(boxmat,1),'r')
% Plot the markers for each box (but I recommend not doing since they are many)
% colors = jet(size(boxmat,1));
% for i = 1:size(boxmat,2)
%     for j = 1:size(boxmat,1)
%         sh = scatter(i, boxmat(j,i),[],colors(j,:),'filled');
%         sh.MarkerFaceAlpha = 0.7;
%         hold on
%     end
% end
xlabel('Time shifts [s]')

