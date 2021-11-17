%% Add paths
clear, clc
addpath('C:/Users/guta_/Documents/MATLAB/eeglab2021.0')
addpath(genpath('./interactiveQRS'))
varsbefore = who; eeglab; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})

runs = {'neurowMIMO_run-01', 'grazMI_run-01', 'neurowMIMO_run-02', 'grazMI_run-02', 'grazME'};

for ru = 1:numel(runs)
    %% Load
    EEG = pop_loadbv('.\', sprintf('sub-03_ses-inside_task-%s.vhdr', runs{ru}), [], []);
    
    %% Correct GA
    Win = 30;
    TR = 1.26;
    thr = 0.02;
    etype = 'Scan Start';    % Name of FMRI slice (slice timing) event.  Unlike fmrib_fastr.m, this is the name of the event.  fmrib_fastr.m takes a vector with slice locations.
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
    
    % Find scans
    EEGMARK = vol_markers_ines(EEG, TR, thr, etype); % sets a Scan Start marker for each volume
    
    % Plot
    figure
    chn = 5;
    scan_lats = [EEGMARK.event(strcmp({EEGMARK.event(:).type},etype)).latency];
    plot(EEGMARK.times(1:scan_lats(1))/1000, EEGMARK.data(chn,1:scan_lats(1)), 'Color', [0.5 0.5 0.5]);
    hold on
    plot(EEGMARK.times(scan_lats(1):scan_lats(end))/1000, EEGMARK.data(chn,scan_lats(1):scan_lats(end)));
    hold on
    plot(EEGMARK.times(scan_lats(end):end)/1000, EEGMARK.data(chn,scan_lats(end):end), 'Color', [0.5 0.5 0.5]);
    hold on
    plot_markers(EEGMARK, {'S  1','S  2','S  5','S  7','S  8','S 10','S 11','S 12'}, [])
    plot_markers(EEGMARK, {etype}, chn)
    xlabel('Time [s]')
    ylabel('Amplitude [uV]')
    title(sprintf('Number of volumes: %d', numel(scan_lats)))
    
    % Trim
    EEGMARK = Trim(EEGMARK, [scan_lats(1), scan_lats(end) + TR * EEGMARK.srate - 1 ], 'idx');
    EEGMARK = Tare(EEGMARK);
    
    % Correct
    [EEGGA, command] = pop_fmrib_fastr(EEGMARK,lpf,L_interp,Win,etype,strig,anc_chk,trig_correct,Volumes,Slices,pre_frac,exc_chan,NPC);
    EEGGA.GAcorrection.win = Win;
    EEGGA.GAcorrection.ANC = anc_chk;
    EEGGA.GAcorrection.interp = L_interp;
    
    % Plot
    figure
    plot(EEGGA.times/1000, EEGGA.data(chn,:));
    hold on
    plot_markers(EEGGA, {'S  1','S  2','S  5','S  7','S  8','S 10','S 11','S 12'}, [])
    plot_markers(EEGGA, {etype}, chn)
    
    %% Downsample
    ds = 250;
    EEGDWS = pop_resample(EEGGA, ds);
    for m = 1:size(EEGDWS.event,2), EEGDWS.event(m).latency = round(EEGDWS.event(m).latency); end
    
    %% Detect QRS
    % Filter
    chn_ecg = 32;
    ECG = pop_select(EEGDWS, 'channel', chn_ecg);
    ECG = pop_eegfiltnew(ECG, 'locutoff',0.5, 'plotfreqz',0);
    ECG = pop_eegfiltnew(ECG, 'hicutoff',30, 'plotfreqz',0);
    EEGDWS.data(chn_ecg, :) = ECG.data(1, :);
    
    % Detect automatically
    freq_band = [4 45]; % [4 45 Hz] to increase QRS detection accuracy (Abreu et al., 201
    reverse = 1;
    [ecg, bpmR, p, R_struct ] = ecgPeakDetection_v4(EEGDWS.data(chn_ecg, :), EEGDWS.srate, freq_band, reverse);
    
    % Plot
    figure
    plot(EEGDWS.times/1000, EEGDWS.data(chn_ecg,:));
    hold on
    plot_markers(EEGDWS, {'S  1','S  2','S  5','S  7','S  8','S 10','S 11','S 12'}, [])
    plot_markers(EEGDWS, {'RT'}, chn)
    
    % Detect interactively
    interactiveQRS(EEGDWS,p);  % Outputs an EEG to the workspace
    input('Continue with interactively QRS-detected EEG? (any key)\n');
    
    %% Correct PA
    % Correct
    EEG.data = double(EEG.data);
    [EEGPA, ~] = pop_fmrib_pas(EEG, 'QRSi', 'obs');
    
    % Plot
    figure
    plot(EEGPA.times/1000, EEGPA.data(chn_ecg,:));
    hold on
    plot_markers(EEGPA, {'S  1','S  2','S  5','S  7','S  8','S 10','S 11','S 12'}, [])
    
    %% Save
    EEG = EEGDWS;
    save(sprintf('sub-03_ses-inside_task-%s_GAcorrected.mat', runs{ru}), 'EEG')
    EEG = EEGPA;
    save(sprintf('sub-03_ses-inside_task-%s_PAcorrected.mat', runs{ru}), 'EEG')
end


