%% Add
clear, clc
addpath('C:/Users/guta_/Documents/MATLAB/eeglab2021.0')
varsbefore = who; eeglab; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})

%% Load (offline-corrected)
% load('sub-03_ses-inside_task-grazMI_run-01_EEGPA.mat','EEG')

%% Set
chan = 'C3';
bands = [1, 4; 4, 8; 8, 13; 13, 30];
band_names = {'Delta', 'Theta', 'Alpha', 'Beta'};
vol_event = 'Vol Start';

%% A) Prep 5 (marks bad TRs)
resamp = 250; 
hp = 1;
lp = 40;
elims = [0, 1.26];
ereject = false;
EEGtr = prep5(EEG,resamp,hp,lp,{vol_event},elims,ereject);
EEG = undo_epochs(EEGtr); % To see bad TRs: EEG.prep_report.('trials_rej_mask')

%% B) Prep 6 (marks bad TRs) (similar to prep 5, skip testing if you wish)
% resamp = 250; 
% hp = 1;
% lp = 40;
% elims = [0, 1.26];
% ereject = false;
% EEGtr = prep6(EEG,resamp,hp,lp,{vol_event},elims,ereject);
% EEG = undo_epochs(EEGtr); % To see bad TRs: EEG.prep_report.('trials_rej_mask')

%% C) Prep 8 (reconstructs bad periods)
resamp = 250;
hp = 1;
lp = 40;
EEG = prep8(EEG,resamp,hp,lp);

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
