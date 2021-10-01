%% Add paths
clear, clc
addpath('C:/Users/guta_/Documents/MATLAB/eeglab2021.0')
addpath(genpath('./interactiveQRS'))
varsbefore = who; eeglab; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})

%% Load
load('../data/sub-03_ses-inside_task-neurowMIMO_run-01.mat','EEG');

%% Correct GA
Win = 30;
TR = 1.26;
thr = 0.02;
etype = 'Vol Start';    % Name of FMRI slice (slice timing) event.  Unlike fmrib_fastr.m, this is the name of the event.  fmrib_fastr.m takes a vector with slice locations.
lpf = [];               % Low pass filter cutoff (default: [ ]=70).
L_interp = 4;           % Interpolation folds (default: [ ]=10).
strig = 1;              % 1 for slice triggers (default) 0 for volume/section triggers
anc_chk = 0;            % 1 to perform ANC 0 No ANC
trig_correct = 0;       % 1 to correct for missing triggers; 0 to NOT correct.
Volumes = [];           % FMRI Volumes.  Needed if trig_correct=1; otherwise use [ ];
Slices = [];            % FMRI Slices/Vol. usage like Volumes.
pre_frac = [];          % Relative location of slice triggers to actual start of slice (slice artifact). default [] = 0.03
exc_chan = 32;
NPC = 'auto';           % Number of principal components to fit to residuals. 0 to skip OBS fitting and subtraction. 'auto' for automatic order selection (default)

[fig, EEGMARK] = vol_markers_ines(EEG, TR, thr, etype); % sets a Scan Start marker for each volume

% Plot
figure
plot_time(EEGMARK, 5, etype)

% GA correct
[EEGGA, command] = pop_fmrib_fastr(EEGMARK,lpf,L_interp,Win,etype,strig,anc_chk,trig_correct,Volumes,Slices,pre_frac,exc_chan,NPC);

EEGGA.GAcorrection.win = Win;
EEGGA.GAcorrection.ANC = anc_chk;
EEGGA.GAcorrection.interp = L_interp;

% Plot
figure
plot_time(EEGGA, 5, etype)

% [EEG,trim1] = GA_vol_based(EEG,c.TR,c.thr,c.lpf,c.L_interp,c.Win,c.etype,c.strig,c.anc_chk,c.trig_correct,c.Volumes,c.Slices,c.pre_frac,c.exc_chan,c.NPC);

%% Downsample
ds = 250;
EEGDWS = pop_resample(EEGGA, ds); % Talvez verificar com função do Matlab
for m = 1:size(EEGDWS.event,2), EEGDWS.event(m).latency = round(EEGDWS.event(m).latency); end

%% Detect QRS
chan_ecg = 32;
EEGECG = pop_select(EEGDWS, 'channel', chan_ecg);
EEGECG = pop_eegfiltnew(EEGECG, 'locutoff',0.5, 'plotfreqz',0);
EEGECG = pop_eegfiltnew(EEGECG, 'hicutoff',30, 'plotfreqz',0);
EEGDWS.data(chan_ecg, :) = EEGECG.data(1, :);

% A - Good detection
ecgchan = 32;
freq_band = [4 45]; % [4 45 Hz] to increase QRS detection accuracy (Abreu et al., 201
reverse = 1;
[ecg, bpmR, p, R_struct ] = ecgPeakDetection_v4(EEGDWS.data(ecgchan, :), EEGDWS.srate, freq_band, reverse);
EEGplot = EEGDWS;
EEGplot.event = R_struct;
%%
plot_time(EEGplot, 32, 'RT')
interactiveQRS(EEGDWS,p);
% or B - Worse detection
% [EEGQRS,~] = pop_fmrib_qrsdetect(EEGDWS,chan_ecg,'QRS','no');
% interactiveQRS(EEGDWS,'QRS');
% or C - No detection
% interactiveQRS(EEGQRS,60);

in = input('Save interactively QRS-detected EEG? (y/n)\n','s');
if strcmp(in,'n'), return, end

%% Correct PA
qrsi_lats = [EEG.event(strcmp({EEG.event(:).type},'QRSi')).latency];
EEG = Trim(EEG, [qrsi_lats(1),qrsi_lats(end)], 'idx');
EEG = Tare(EEG);
EEG.data = double(EEG.data);

[EEGPA, ~] = pop_fmrib_pas(EEG, 'QRSi', 'obs');

figure
plot_time(EEGPA, 5, [])




