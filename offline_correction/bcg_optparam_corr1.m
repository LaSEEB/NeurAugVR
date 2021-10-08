addpath('/home/iesteves/MR_Correction_Pipeline/eeglab2019_0') 
eeglab

%%  PA Correction using optimal parameters

cd '/strombolihome/iesteves/MR_Correction_Pipeline/Preprocessing/ICA'
files = dir('*resting*.set');

filenames = {files.name};
filefolder = {files.folder};

param_path = '/strombolihome/iesteves/MR_Correction_Pipeline/Preprocessing/';
path = '/strombolihome/iesteves/MR_Correction_Pipeline/Preprocessing/';

locfilepath = '/strombolihome/iesteves/MR_Correction_Pipeline/Pipeline';
locfilename = strcat(locfilepath, '/chan_loc.ced'); %'Chan_loc_x20t.ced';

chan = 6; % change channel number for plots
ecgchan = 32;

%% Load EEG file

wbkg_all = [0, 0.1, 0.3, 0.5, 0.7, 1]; % physiological background weight: varies from 0 to 1, in steps of 0.1 
%wbkg_all = [0, 0.5, 1]; % physiological background weight: varies from 0 to 1, in steps of 0.1 

for n = 3:length(filenames)
    
    % subject and task
    strfile = strsplit(filenames{n}, {'_', '.'});
    subject = strfile{2};
    task = strfile{1};

    % load file 
    [EEGICA, ~] = pop_loadset(filenames{n}, [filefolder{n}, '/']);

    % load channel location
    EEGICA.chanlocs = readlocs(locfilename);

    % load optimal parameters
    load([param_path, 'opt_param_restIAB_', task, '_', subject, '.mat'])
    %load([param_path, 'opt_param_', task, '_', subject, '.mat'])

    %% PA correction Rodolfo Abreu plugin
    %% 6 - Applying EEG BCG Correction

    fprintf('\n--->>>---> 6 - Applying BCG Correction...\n\n');
    % Reference: https://github.com/rmabreu/BCG_Artefact_Correction

    % change the following parameters accordingly
    art_harm = 5;           % number of harmonics to include in the BCG artefact correction quantification
    win_hz = 0.065;         % window length for which the BCG artefact correction will be assessed
    filters = [ 0.5 45 ];   % low and high cutoff values of the band-pass filtering of EEG data
    TR = 1.26;           % repetition time of the fMRI acquisition [s]
    harm_thr = 0.25;        % threshold for inclusion of harmonics {recommended = 0.25; for higher variability across subjects and channels, the threshold should be lower}
    plt = 0;                % 1/0 = do/don't display the results
    nb_plts = 4;            % number of subplots within each plot for displaying the selected IC time-courses


    % load all necessary matrices (data, activations and mixing matrix from ICA)

    dataset = EEGICA;
    sph = EEGICA.icasphere;
    wghts = EEGICA.icaweights;
    activations = EEGICA.icaact;

    mix_matrix = wghts * sph;

    Kp = EEGICA.qrs.allpeaks;


    dataset.Kp = unique(sort(round(Kp)));   
    dataset.ecg = dataset.data(32,:);  
    dataset.data = dataset.data(1:31,:);

    % retrieve R peak annotations and compute appropriate limits for epoching
    % EEG data using the R peaks as triggers (dependent on the subjects' heart
    % rate)
    R = dataset.Kp; 
    m_R = (min(diff(R)) / dataset.srate) * 1000;
    %lim_inf_ms = -100;
    lim_inf_ms = -50;
    if m_R < 4 * abs(lim_inf_ms)
        lim_sup_ms = lim_inf_ms + (4 * abs(lim_inf_ms));
    else
        lim_sup_ms = lim_inf_ms + m_R;
    end
    limits_ecg = [ lim_inf_ms lim_sup_ms ];

    % define the standard time-delay (210 ms) between R peak and BCG artefact
    % occurrences (according to Allen et al., NeuroImage 1998)
    delay_QRS = 0.21;

    for wbkg = wbkg_all    
        % PROJIC-AAS
        % n_win = numbers of averaging windows
        % k_clusters = number of clusters for ECG-related ICs
       
        % k = opt_param_IAB(cmp(opt_param_IAB.weight, wbkg),:).OBS_k; % optimal k for weight w
        % win = opt_param_IAB(cmp(opt_param_IAB.weight, wbkg),:).OBS_pc;  % optimal win for weight w

        k = opt_param_IAB(cmp(opt_param_IAB.weight, wbkg),:).OBS_k; % optimal k for weight w
        win = opt_param_IAB(cmp(opt_param_IAB.weight, wbkg),:).OBS_pc;  % optimal win for weight w


        k_clusters = k; n_win = win; 
        [ ~, ~, eeg_bcg ] = BCG_Correction_PROJIC_OBS(dataset, activations, ...
           mix_matrix, limits_ecg, filters, n_win, k_clusters, art_harm, harm_thr, ...
            TR, win_hz, delay_QRS, plt, nb_plts);


          %k = opt_param(opt_param.weight == wbkg,:).OBS_k; % optimal k for weight w
          %win = opt_param(opt_param.weight == wbkg,:).OBS_pc;  % optimal win for weight w

    %     k_clusters = k; npc = win; 
    %     [ ~, ~, eeg_bcg ] = BCG_Correction_PROJIC_OBS(dataset, ...
    %     activations, mix_matrix, limits_ecg, filters, npc, k_clusters, art_harm, ...
    %     harm_thr, TR, win_hz, delay_QRS, plt, nb_plts);


        EEGPAR = EEGICA;
        EEGPAR.data(1:end-1, :) = eeg_bcg; % eeg backreconstructed without artifact
        EEGPAR.PAcorrection.method = 'projic_obs';
        EEGPAR.PAcorrection.weight = wbkg;
        EEGPAR.PAcorrection.k = k_clusters;
        EEGPAR.PAcorrection.win = n_win;

        % plot
%         N = length(EEGPAR.data(chan, :));
% 
%         fft_qrs = (abs(fft(EEGICA.data(chan,:))).^2)/N;
%         fft_par = (abs(fft(EEGPAR.data(chan,:))).^2)/N;
%         fft_ecg = (abs(fft(EEGICA.data(ecgchan,:))).^2)/N;
%         L_qrs = length(fft_qrs);
%         ymax = max(fft_qrs);
%         fvalue = EEGICA.srate*(0:(L_qrs/2))/L_qrs;
%         [~, f_40] = min(abs(fvalue-40)); % find the index that corresponds to 40 Hz
%         [~, f_01] = min(abs(fvalue-0.10)); % find the index that corresponds to 0.1 Hz
% 
%         figure 
%         subplot(2,1,1);
%         plot((1:length(EEGICA.data(chan,:)))/EEGICA.srate, EEGICA.data(ecgchan,:))
%         hold on
%         plot((1:length(EEGICA.data(chan,:)))/EEGICA.srate, EEGICA.data(chan,:))
%         plot((1:length(EEGPAR.data(chan,:)))/EEGPAR.srate, EEGPAR.data(chan,:))
% 
%         titleplot = strcat({['Pulse Artifact Correction', ' - ', subject],['Channel: ', EEGPAR.chanlocs(chan).labels, ...
%             '; wbkg=', num2str(wbkg)]});
%         title(titleplot);
%         ylabel(' Amplitude[uV] ')
%         xlabel('Time[s]')
%         legend('ECG', 'Original EEG', 'PA Removed EEG');
%         hold off
% 
%         subplot(2,1,2); 
%         plot(fvalue(f_01:f_40), fft_ecg(f_01:f_40))
%         hold on
%         plot(fvalue(f_01:f_40), fft_qrs(f_01:f_40))
%         plot(fvalue(f_01:f_40), fft_par(f_01:f_40))
%         title('FFT PA Removed EEG')
%         ylabel('Spectrum Amplitude')
%         xlabel('Frequency [Hz]')
%         legend('ECG', 'Original EEG', 'PA Removed EEG');
%         hold off
%         print([path, '/Figures/', task, '_', subject, '_PA_filtered_OBS_wbkg_', num2str(wbkg)], '-dpng')
%         savefig([path, '/Figures/', task, '_', subject, '_PAfiltered_OBS_wbkg_', num2str(wbkg),'.fig'])
% 
%         pop_eegplot(EEGPAR, 1, 1, 1);
%         figure; pop_spectopo(EEGPAR, 1, [0 1000*EEGPAR.pnts/EEGPAR.srate], 'EEG' , 'freq', [6 16 32], 'freqrange',[1 40],'electrodes','off');
%         savefig([path, '/Figures/', task, '_', subject, '_PAfiltered-projicOBS_spectopo_wbkg_', num2str(wbkg), '.fig'])

%        PAfilename = [task, '_', subject, '_PAfiltered-IABprojicOBS_wbkg_', num2str(wbkg), '.set'];
        PAfilename = [task, '_', subject, '_PAfiltered-restIABprojicOBS_wbkg_', num2str(wbkg), '.set'];
        PAfilepath = [path, '/PAfiltered_PROJIC_OBS/'];
        pop_saveset(EEGPAR, 'filename', PAfilename, 'filepath',PAfilepath);

    end
end