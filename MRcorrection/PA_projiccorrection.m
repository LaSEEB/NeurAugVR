
function EEGPAR = PA_projiccorrection(EEGICA, PAmethod, opt_param, wbkg, TR)

    % change the following parameters accordingly
    art_harm = 10;           % number of harmonics to include in the BCG artefact correction quantification
    win_hz = 0.065;         % window length for which the BCG artefact correction will be assessed
    filters = [0.5 40];   % low and high cutoff values of the band-pass filtering of EEG data
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
    lim_inf_ms = -100;
    if m_R < 4 * abs(lim_inf_ms)
        lim_sup_ms = lim_inf_ms + (4 * abs(lim_inf_ms));
    else
        lim_sup_ms = lim_inf_ms + m_R;
    end
    limits_ecg = [ lim_inf_ms lim_sup_ms ];

    % define the standard time-delay (210 ms) between R peak and BCG artefact
    % occurrences (according to Allen et al., NeuroImage 1998)
    delay_QRS = 0.21;

    switch PAmethod
       case 'OBS'
            k_clusters = opt_param(cmp(opt_param.weight, wbkg),:).OBS_k; % optimal k for weight w
            npc = opt_param(cmp(opt_param.weight, wbkg),:).OBS_pc;  % optimal win for weight w

            [ ~, ~, eeg_bcg ] = BCG_Correction_PROJIC_OBS(dataset, activations, ...
               mix_matrix, limits_ecg, filters, npc, k_clusters, art_harm, harm_thr, ...
                TR, win_hz, delay_QRS, plt, nb_plts);

            EEGPAR = EEGICA;
            EEGPAR.data(1:end-1, :) = eeg_bcg; % eeg backreconstructed without artifact
            EEGPAR.PAcorrection.method = 'projic_obs';
            EEGPAR.PAcorrection.weight = wbkg;
            EEGPAR.PAcorrection.k = k_clusters;
            EEGPAR.PAcorrection.npc = npc;
        
       case 'AAS'
            k_clusters = opt_param(cmp(opt_param.weight, wbkg),:).AAS_k; % optimal k for weight w
            n_win = opt_param(cmp(opt_param.weight, wbkg),:).AAS_win;  % optimal win for weight w

            [ ~, ~, eeg_bcg ] = BCG_Correction_PROJIC_AAS(dataset, activations, ...
               mix_matrix, limits_ecg, filters, n_win, k_clusters, art_harm, harm_thr, ...
                TR, win_hz, delay_QRS, plt, nb_plts);

            EEGPAR = EEGICA;
            EEGPAR.data(1:end-1, :) = eeg_bcg; % eeg backreconstructed without artifact
            EEGPAR.PAcorrection.method = 'projic_aas';
            EEGPAR.PAcorrection.weight = wbkg;
            EEGPAR.PAcorrection.k = k_clusters;
            EEGPAR.PAcorrection.win = n_win;
    end
