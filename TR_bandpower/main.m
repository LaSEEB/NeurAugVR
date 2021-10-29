%% Add
clear, clc
addpath('C:/Users/guta_/Documents/MATLAB/eeglab2021.0')
varsbefore = who; eeglab; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})

%% Load (simple-offline-corrected + preprocessed-4)
load('sub-3_ses-2_run-1_EEGpre4_wreport.mat', 'EEG');

%% Set
chan = 'C3';
bands = [1, 4; 4, 8; 8, 13; 13, 30];
band_names = {'Delta', 'Theta', 'Alpha', 'Beta'};
vol_event = 'Vol Start';

%% A) Use bandpower() on each TR
% Calculate
[tr_vecs, pow_vecs] = tbandpower(EEG, chan, vol_event, bands);

%% B) Bandpass signal and then calculate power for each TR
% Calculate
clc
[tr_vecs, pow_vecs] = tfiltpower(EEG, chan, vol_event, bands);

%% Plot
figure
for b = 1:size(bands,1)
    subplot(numel(band_names), 1, b)
    plot(tr_vecs(1,:), pow_vecs(b,:))
    ylabel(band_names{b})
    if b == 1, title('TR-bandpower'); end
end
xlabel('TR')
