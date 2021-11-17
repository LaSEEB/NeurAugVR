%% Add
clear, clc
addpath('C:/Users/guta_/Documents/MATLAB/eeglab2021.0')
varsbefore = who; eeglab; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})

%% Load (offline-corrected)
% load('sub-03_ses-inside_task-grazMI_run-01_EEGPA.mat','EEG')
load('sub-04_ses-inside_task-grazMIMO_run-0001_PAcorrected-obs.mat','EEGPAR') % fix 1: EEG -> EEGPAR
EEG = EEGPAR;

%% Set
chan = 'C3';
bands = [1, 4; 4, 8; 8, 13; 13, 30];
band_names = {'Delta', 'Theta', 'Alpha', 'Beta'};
vol_event = 'Scan Start'; % 'Vol Start'/'Scan Start'

%% A) Prep 5 (marks bad TRs)
resamp = 250; 
hp = 1;
lp = 40;
elims = [0, 1.26 - 1/EEG.srate];
ereject = false;
EEGtr = prep5(EEG,resamp,hp,lp,{vol_event},elims,ereject);
EEGtr = remove_event_repetitions(EEGtr); % If epochs overlap and events get counted twice (should not happen if epoching interval is correct), remove repetitions
EEG = undo_epochs(EEGtr); % To see bad TRs: EEG.prep_report.('trials_rej_mask')

%% B) Prep 6 (marks bad TRs) (similar to prep 5, skip testing if you wish)
% resamp = 250; 
% hp = 1;
% lp = 40;
% elims = [0, 1.26 - 1/EEG.srate];
% ereject = false;
% EEGtr = prep6(EEG,resamp,hp,lp,{vol_event},elims,ereject);
% EEGtr = remove_event_repetitions(EEGtr); % If epochs overlap and events get counted twice (should not happen if epoching interval is correct), remove repetitions
% EEG = undo_epochs(EEGtr); % To see bad TRs: EEG.prep_report.('trials_rej_mask')

%% C) Prep 8 (reconstructs bad periods)
resamp = 250;
hp = 1;
lp = 40;
EEG = prep8(EEG,resamp,hp,lp);

%% A) Use bandpower() on each TR
% Calculate
[tr_vecs, pow_vecs] = tbandpower_tr(EEG, chan, vol_event, bands);

%% B) Bandpass signal and then calculate power for each TR
% Calculate
[tr_vecs, pow_vecs] = tfiltpower_tr(EEG, chan, vol_event, bands);

%% C) Band-task mean ERSP
% Slide through TR's
% onsets = [EEG.event(strcmp({EEG.event(:).type},'Scan Start')).latency];
% trial = [-4, 4] * 1.26;  % s

% Slide through every point
onsets = 1:size(EEG.data,2);
trial = [-5, 5];  % s

time_id = 1.75;  % Each ERSP will be saved in the segment at time_id s (e.g.: 1.75s) inside the trial (0s represents the cross appearing)

[x_vecs, ersp_vecs] = tfERSP(EEG, chan, onsets, bands, trial, time_id);

%% Plot
plot_vecs = ersp_vecs;
figure
for b = 1:size(bands,1)
    subplot(numel(band_names), 1, b)
    plot(x_vecs(1,:), plot_vecs(b,:))
    ylabel(band_names{b})
    if b == 1, title('Bandpowers'); end
end
xlabel('Onset id')
