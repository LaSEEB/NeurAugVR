%% Add paths
clear, clc
addpath('C:/Users/guta_/Documents/MATLAB/eeglab2021.0')
varsbefore = who; eeglab; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})

%% Load (simple-offline-corrected + preprocessed-4)
load('EEGnocrop_pre4_wreport.mat', 'EEG');
% load('EEGpre4_wreport.mat', 'EEG');
%% Set channel, bands and time intervals (indexes)
chan = 'C3';
chn = find(ismember({EEG.chanlocs(:).labels}, chan));
bands = [1, 4; 4, 8; 8, 13; 13, 30];
band_names = {'Delta', 'Theta', 'Alpha', 'Beta'};
vol_lats = [EEG.event(strcmp({EEG.event(:).type},'Vol Start')).latency];

%% A) Use bandpower() to calculate bandpower for each TR
pow_vecs = zeros(size(bands, 1),size(vol_lats,2)-1);
tr_vecs = 1:numel(vol_lats)-1;
tr_vecs(2,:) = (EEG.times(vol_lats(tr_vecs)) - EEG.times(vol_lats(1)))/1000;
figure
for b = 1:size(bands,1)
    for v = 1:numel(vol_lats)-1
        p = bandpower(EEG.data(chn,vol_lats(v):vol_lats(v+1)),EEG.srate,bands(b,:));
        pow_vecs(b,v) = p;
    end
    pow_vecs(b,:) = -1 + 2.*(pow_vecs(b,:) - min(pow_vecs(b,:)))./(max(pow_vecs(b,:)) - min(pow_vecs(b,:)));
    subplot(numel(band_names), 1, b)
    plot(tr_vecs(1,:), pow_vecs(b,:))
    ylabel(band_names{b})
end
xlabel('TR')
subplot(numel(band_names), 1, 1)
title('time-bandpower')
% save('sub-03_ses-inside_task-neurowMIMO_run-01_tbandpower.mat', 'tr_vecs', 'pow_vecs')

%% B) Use EEGLAB newtimef() to calculate bandpower for each TR
baseline = nan; % Indiferent actually
scale = 'abs';
cycles = 0;
ylims = [0, 30];

figure
nv = 0;  % Margin TR's (e.g., nv = 2 => 5 TR's)
timelims = [EEG.times(vol_lats(1)) EEG.times(vol_lats(1+nv*2+1))] - EEG.times(vol_lats(1));
tr_vecs = 1+nv : numel(vol_lats)-nv-1;
tr_vecs(2,:) = (EEG.times(vol_lats(tr_vecs)) - EEG.times(vol_lats(1)))/1000;
amp_vecs = zeros(size(bands, 1),size(tr_vecs,2));

for i = 1:size(tr_vecs,2)
    v = tr_vecs(1, i);
    data_vec = EEG.data(chn,vol_lats(v-nv):vol_lats(v+nv+1));
    [~, ~, ~, time_vec, freq_vec, ~, ~, amp_mat] = newtimef_trueamp(data_vec, size(data_vec,2), timelims, EEG.srate, cycles, 'freqs', ylims,'plotitc','off', 'baseline', baseline, 'scale', scale, 'plotersp', 'off');
    amp_mat = double(amp_mat);
    
    for b = 1:size(bands,1)
        freqmask = freq_vec >= bands(b,1) & freq_vec <= bands(b,2);
        amp_vecs(b,i) = mean(amp_mat(freqmask,:),'all');
    end
end
for b = 1:size(bands,1)
    amp_vecs(b,:) = -1 + 2.*(amp_vecs(b,:) - min(amp_vecs(b,:)))./(max(amp_vecs(b,:)) - min(amp_vecs(b,:)));
    subplot(numel(band_names), 1, b)
    plot(tr_vecs(1, :), amp_vecs(b,:))
    ylabel(band_names{b})
end
xlabel('TR')
subplot(numel(band_names), 1, 1)
title('tf-bandpower-eachTR')
% save('sub-03_ses-inside_task-neurowMIMO_run-01_tfbandpower-eachTR.mat', 'tr_vecs', 'amp_vecs')

