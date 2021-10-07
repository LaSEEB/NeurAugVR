eeglabfolder1 = '/home/mfleury/POSTDOC/LIBRAIRY/eeglab2019_0';
eeglabfolder2 = '/home/mfleury/POSTDOC/LIBRAIRY/EEGLAB';

addpath(eeglabfolder1) 
%eeglab 

%%

[~, name] = system('hostname');
name = lower(name);

if strcmp(name(1:end-1), 'stromboli')
    path = '/home/iesteves/DataMR/Preprocessing';
    files = dir('/home/iesteves/DataMR/RawData/sub*020/*.vhdr');
    chan_loc = '/home/iesteves/MR_Correction_Pipeline/Pipelmatine/chan_loc.ced';
elseif strcmp(name(1:end-1), 'sesimbra')
    path = '/home/mfleury/POSTDOC/DATA/STUDY_NeurAugVr/pre-processing/gr_sujetssains/su_01/eeg/';
    files = dir('/home/mfleury/POSTDOC/DATA/STUDY_NeurAugVr/ses-inside/su_01/eeg/*.vhdr');
    chan_loc = '/home/mfleury/POSTDOC/LIBRAIRY/MR_Correction_Pipeline/Pipeline/chan_loc.ced';
end

cd '/home/mfleury/POSTDOC/DATA/STUDY_NeurAugVr' % CHANGE ACCORDINGLY

filenames = {files.name};
filefolder = {files.folder};

locfilepath = '/home/mfleury/POSTDOC/LIBRAIRY/MR_Correction_Pipeline/Pipeline';
locfilename = strcat(locfilepath, '/chan_loc.ced'); %'Chan_loc_x20t.ced';

chan = 6; % change channel number for plots


TR = 1.26; % repetition time in seconds
thr = 0.02; % default threshold to find volume markers

checkfolders(path);
%% Load EEG file

% Convert from BrainVision original files to .mat files, without any
% processing
for n = 1:length(filenames)
    EEG = pop_loadbv([filefolder{n}, '/'], filenames{n});
    
    filename = filenames{n}(1:end-5);
    EEG.setname = filename;
    save([path,'/Original_mat/',filename], 'EEG');

end

%%
for n = 1:length(filenames)
    

    % load file 
    filename = filenames{n}(1:end-5);
    EEGfile = load([path,'/Original_mat/',filename]);

    EEG = EEGfile.EEG;
    
    % load channel location
    EEG.chanlocs = readlocs(locfilename);

 
    % Read acquisition details: subject, session, task, run
    [subject, session, task, run] = acquisition_names(filename);
    if any(strcmp(task, {'task-eegcalibrationout', 'task-breathholdout'}))
        continue
    end

    % change event types from number to string
    event_types_char = cellfun(@num2str, {EEG.event(:).type}, 'UniformOutput', false);
    [EEG.event(:).type] = deal(event_types_char{:});
 
    %% Scanner markers
    disp('--->>>---> 1 - Placing volume markers...');

    [fig, EEGMARK, vols] = volume_markers(EEG, TR, thr); % sets a Scan Start marker for each volume 
    EEGMARK.volmarkers = vols;

    save([path, '/Markers/', filename, '_markers'], 'EEGMARK')

    %% GA correction
    disp('--->>>---> 2 - Removing the gradient artifact...');

    %EEGMARK = pop_eegfiltnew(EEGMARK, 1, []);

    EEGMARK.data = double(EEGMARK.data);

    % Apply artifact correction

    lpf = 70;               % Low pass filter cutoff (default: [ ]=70).
    L_interp = 4;                  % Interpolation folds (default: [ ]=10).
    Win = 30;               %  Number of artifacts in avg. window (default: [ ]=30)
    etype = 'Scan Start';   % Name of FMRI slice (slice timing) event.  Unlike
                            % fmrib_fastr.m, this is the name of the event.  fmrib_fastr.m 
                            % takes a vector with slice locationsInvalid MEX-file '/strombolihome/iesteves/eeglab2019_0/plugins/fmrib1.21/fastranc.mexa64'.
    strig = 1;              % 1 for slice triggers (default)
                            % 0 for volume/section triggers
    anc_chk = 0;            % 1 to perform ANC
                            % 0 No ANC
    trig_correct = 0;       % 1 to correct for missing triggers;
                            % 0 to NOT correct.
    Volumes = [];           % FMRI Volumes.  Needed if trig_correct=1; 
                            % otherwise use [ ];
    Slices = [];            % FMRI Slices/Vol. usage like Volumes.
    pre_frac = 0;          % Relative location of slice triggers to actual start
                            %of slice (slice artifact). default [] = 0.03
    exc_chan = 32;          
    NPC = 'auto';           % Number of principal components to fit to residuals.
                            % 0 to skip OBS fitting and subtraction.
                            % 'auto' for automatic order selection (default)

    [EEGGAR, command] = pop_fmrib_fastr(EEGMARK,lpf,L_interp,Win,etype,...
        strig,anc_chk,trig_correct,Volumes,Slices,pre_frac,exc_chan,NPC); 

    EEGGAR.GAcorrection.win = Win;
    EEGGAR.GAcorrection.ANC = anc_chk;
    EEGGAR.GAcorrection.interp = L_interp;

    event_types = {EEGGAR.event(:).type};
    vol_index = find(strcmp('Scan Start', event_types));
    first_vol = EEGGAR.event(vol_index(1)).latency;
    last_vol = EEGGAR.event(vol_index(end)).latency;


    %% Delete signal portion outside the correction (which corresponds to fMRI acquisition)
    disp('--->>>---> 3 - Cutting the signal to match fMRI acquisition...');

    EEGCUT = pop_select(EEGGAR, 'point', [first_vol last_vol+TR*EEGGAR.srate-1]);
    
    for ev = 1:length(EEG.event)
        EEGCUT.event(ev).latency = round(EEGCUT.event(ev).latency);
    end
    
    %% Downsampling the signal from 5000Hz to 250Hz
    disp('--->>>---> 4 - Downsampling the signal to 250 Hz...');

    newfreq = 250;

    EEGDWS = pop_resample(EEGCUT, newfreq);

    save([path, '/DWS/', filename, '_dws'], 'EEGDWS')

    %% QRS detection
    disp('--->>>---> 5 - QRS detection: STEP 1 - automatic R-peak detection ...');

    ecgchan = 32;
    fs = EEGDWS.srate;
    signal = EEGDWS;
    freq_band = [4 45]; % [4 45 Hz] to increase QRS detection accuracy (Abreu et al., 201
    reverse = 1;
    [ecg, bpmR, p, R_struct ] = ecgPeakDetection_v4(signal.data(ecgchan, :), fs, freq_band, reverse);

    mkdir([path,'/qrs_points'])
    save([path,'/qrs_points/qrspoints_',filename,'.mat'], 'p')

end

%%
for n = 1:length(filenames)
    
    filename = filenames{n}(1:end-5);
    
    % Read acquisition details: subject, session, task, run
    [subject, session, task, run] = acquisition_names(filename);
    if any(strcmp(task, {'task-eegcalibrationout', 'task-breathholdout'}))
        continue
    end
     %% QRS detection
    disp('--->>>---> 5 - QRS detection: STEP 2 - manual R-peak verification ...');
    
    EEGDWSfile = load([path, '/DWS/', filename, '_dws']);
    EEGDWS = EEGDWSfile.EEGDWS;
    
    qrspoints = load([path,'/qrs_points/qrspoints_',filename,'.mat']);
    starter_marker_lats = qrspoints.p;
    
    interactiveQRS(EEGDWS, starter_marker_lats);
    pause
    EEGQRS = EEG;
    savefig([path, '/Figures/rpeaks_',filename]);
    
    event_types = {EEGQRS.event(:).type};
    qrs_ev = find(strcmp('QRSi', event_types));
    qrs_index = round([EEGQRS.event(qrs_ev).latency]);
    bpm = round((60*length(qrs_index))/(EEGQRS.pnts/EEGQRS.srate));

    EEGQRS.qrs.allbpm = bpm;
    EEGQRS.qrs.allpeaks = qrs_index;
    EEGQRS.qrs.allpeaks_nr = length(qrs_index);

    save([path, '/QRS/', filename, '_qrs'], 'EEGQRS')
end

%% ICA
ecg_chan = 32;
for n = 1:length(filenames)
    disp('--->>>---> 6 - ICA weights for PA PROJIC correction ...');
    filename = filenames{n}(1:end-5);
    % Read acquisition details: subject, session, task, run
    [subject, session, task, run] = acquisition_names(filename);
    if any(strcmp(task, {'task-eegcalibrationout', 'task-breathholdout'}))
        continue
    end
    
    EEGQRSfile = load([path, '/QRS/', filename, '_qrs']);
    EEGQRS = EEGQRSfile.EEGQRS;
    
    pop_editoptions('option_computeica', 1 );
    EEGICA = pop_runica(EEGQRS, 'icatype', 'runica', 'chanind', 1:ecg_chan, 'extended', 1);

    save([path, '/ICA/', filename, '_ica'], 'EEGICA')
end

%% PA correction
powerfreqfolder = [path, '/PowerFreq'];

for n = 1:length(filenames)
    
    disp('--->>>---> 7 - Applying PA Correction...');
    filename = filenames{n}(1:end-5);
        % Read acquisition details: subject, session, task, run
    [subject, session, task, run] = acquisition_names(filename);
    if any(strcmp(task, {'task-eegcalibrationout', 'task-breathholdout'}))
        continue
    end
    
    load([path, '/ICA/', filename, '_ica']);
    
    % Power in the artifact/physiological background windows for PROJIC,
    % PROJIC-OBS and PROJIC-AAS using diferent parameters
    projic_k = 2:12;
    obs_k = 2:12;
    obs_npc = 3:12;
    aas_k = 2:12;
    aas_win = 10:10:50;
    PA_powerfreq(EEGICA, filename, powerfreqfolder, projic_k, obs_k, obs_npc, aas_k, aas_win, TR);

    % C ratio computed for each physiological weight and parameters with
    % maximum C ratio for each weight
    plt = 0;
    opt_param = PA_Cratio(powerfreqfolder, filename, projic_k, obs_k, obs_npc, aas_k, aas_win, plt);
    
    % PA correction using the optimal parameters for the chosen weight and
    % method
    wbkg = 1;
    PAmethod = 'OBS';
    EEGPAfiltered = PA_projiccorrection(EEGICA, PAmethod, opt_param, wbkg, TR);
    save([path, '/PAfiltered/', filename, '_pafiltered'], 'EEGPAfiltered')
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-processing

%% EEG filter 0.3-70Hz + bireref + ICA weights
% nchan  = 31;
% lp = 0.3;
% hp = 70;
% ecgchan = 32;
% for n = 1:length(filenames)
%     
%     disp('--->>>---> 8 - Applying Filter + Re-referencing + Computing ICA weights...');
%     filename = filenames{n}(1:end-5);
%     
%     % Read acquisition details: subject, session, task, run
%     [subject, session, task, run] = acquisition_names(filename);
%     if any(strcmp(task, {'task-eegcalibrationout', 'task-breathholdout'}))
%         continue
%     end
%     
%     load([path, '/PAfiltered/', filename, '_pafiltered']);
%     
%     % Filtering
%     EEGFILT = filter_eeglab(EEGPAfiltered, lp, hp);
%     EEGFILT.bpfilter = [lp  hp];
%     
%     % Re-referencing
%     EEGREF = EEGFILT;
%     
%     el = EEGREF.data([1:ecgchan-1], :);
%     [ bi_mean, ~] = myBiweight(el');
% 
%     EEGREF.data([1:ecgchan-1], :) = el - repmat(bi_mean, nchan, 1); % re-referencing
%     EEGREF.bimean = bi_mean;
%     
%     % Computing IC weights
%     EEGREF = pop_runica(EEGREF, 'icatype', 'runica', 'chanind', [1:ecgchan-1], 'extended', 1);
%     save([path, '/Reref/', filename, '_reref'], 'EEGREF')
% 
% end
%%

%% IClabel IC rejection
% rmpath(eeglabfolder1)
% addpath(eeglabfolder2)
% eeglab
% 
% nchan  = 31;
% channels = 1:31;
% for k = 1:length(filenames)
%     
%     disp('--->>>---> 9 - ICA denoising + Outlier removal...');
%     filename = filenames{k}(1:end-5);
%     
%     % Read acquisition details: subject, session, task, run
%     [subject, session, task, run] = acquisition_names(filename);
%     if any(strcmp(task, {'task-eegcalibrationout', 'task-breathholdout'}))
%         continue
%     end
%     
%     load([path, '/Reref/', filename, '_reref'])
%     
%     EEG = EEGREF;
%     EEG.data = double(EEG.data);
%     % icalabel
%     EEG = iclabel(EEG);
% 
%     % classifications
%     EEG.etc.ic_classification.ICLabel.classifications % n*7 matrix where n=number of IC's and 7=number o classes
% 
%     % show classes
%     EEG.etc.ic_classification.ICLabel.classes
%     % ans = 'Brain'    'Muscle'    'Eye'    'Heart'    'Line Noise'    'Channel Noise'    'Other'
% 
%     %IC threshold value in %
%     thres = 0.9; % 50%
%     
%     % Find IC's to reject
%     for n = 1:length(EEG.etc.ic_classification.ICLabel.classifications)
% 
%         % find the index of the maximum classification score for each IC
%         [val, idx] = max(EEG.etc.ic_classification.ICLabel.classifications(n, :));
%         % if val > 50% and idx not 1 (Brain), then add to list.
%         % todo:check also brain percentage
%         if (val > thres) && (all(idx ~= [1, 7]))
%             % add to a rejection list
%             rej_ic(:,n) = n; % add the IC number
%         else
%             rej_ic(:,n) = NaN;
%         end
% 
%     end
% 
%     rej_ic = rej_ic(~isnan(rej_ic)); %clear the NaNs
%     
%     % reject components
%     EEG = pop_subcomp(EEG, rej_ic, 0);
%     EEG.rej_ic = rej_ic;
%     save([path, '/ICremoved/', filename, '_icremoved'], 'EEG')
%     
%     signal = EEG.data(channels, :);
%     [newsig_capmean, ~, ~] = outlier_rejection(signal, 'capmean');
% 
%     EEGCAP = EEG; EEGCAP.data(channels, :) = newsig_capmean;
%     EEGCAP.outlier.method = 'cap';
%     EEGCAP.outlier.win = '1';
%     EEG = EEGCAP;
%     save([path, '/Preproc/', filename, '_preproc'], 'EEG')

%end
