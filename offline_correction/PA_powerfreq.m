function PA_powerfreq(EEGICA, filename, powerfreqfolder, projic_k, obs_k, obs_npc, aas_k, aas_win, TR)
% Reference: https://github.com/rmabreu/BCG_Artefact_Correction
   
    % change the following parameters accordingly
    art_harm = 10;           % number of harmonics to include in the BCG artefact correction quantification
    win_hz = 0.065;         % window length for which the BCG artefact correction will be assessed
    filters = [ 0.5 40];   % low and high cutoff values of the band-pass filtering of EEG data
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
    lim_inf_ms = -100; % -100, altered to -50 following Rodolfo's suggestion, since 
        % we are detecting the Q of the QRS instead of the R-peak;
    if m_R < 4 * abs(lim_inf_ms)
        lim_sup_ms = lim_inf_ms + (4 * abs(lim_inf_ms));
    else
        lim_sup_ms = lim_inf_ms + m_R;
    end
    limits_ecg = [ lim_inf_ms lim_sup_ms ];

    % define the standard time-delay (210 ms) between R peak and BCG artefact
    % occurrences (according to Allen et al., NeuroImage 1998)
    delay_QRS = 0.21;

    % PROJIC
    k_clusters = projic_k; % 2:12; % k_clusters = number of clusters
    [ PowerFreq, PowerFreq_bkg ] = BCG_Correction_PROJIC(dataset, activations, ...
        mix_matrix, limits_ecg, filters, k_clusters, art_harm, harm_thr, TR, win_hz, ...
        plt, nb_plts);
    save([powerfreqfolder, '/PowerFreq','_',filename,'.mat'], 'PowerFreq')
    save([powerfreqfolder, '/PowerFreq_bkg','_', filename,'.mat'], 'PowerFreq_bkg')

    %PROJIC-OBS
    k_clusters = obs_k; % 2:12; % vector with number of clusters;
    npc = obs_npc; % 3:12;  % vector with number of principal components
    [ PowerFreq_obs, PowerFreq_bkg_obs ] = BCG_Correction_PROJIC_OBS(dataset, activations, ...
        mix_matrix, limits_ecg, filters, npc, k_clusters, art_harm, harm_thr, ...
        TR, win_hz, delay_QRS, plt, nb_plts);
    save([powerfreqfolder, '/PowerFreq_obs','_', filename,'.mat'], 'PowerFreq_obs')
    save([powerfreqfolder, '/PowerFreq_bkg_obs','_', filename,'.mat'], 'PowerFreq_bkg_obs')


    % PROJIC-AAS
    k_clusters = aas_k; % 2:12; % vector with number of clusters;
    n_win = aas_win; % 10:10:50; % vector with number of averaging windows;
    [ PowerFreq_aas, PowerFreq_bkg_aas ] = BCG_Correction_PROJIC_AAS(dataset, activations, ...
        mix_matrix, limits_ecg, filters, n_win, k_clusters, art_harm, harm_thr, ...
        TR, win_hz, delay_QRS, plt, nb_plts);
    save([powerfreqfolder, '/PowerFreq_aas','_', filename,'.mat'], 'PowerFreq_aas')
    save([powerfreqfolder, '/PowerFreq_bkg_aas','_', filename,'.mat'], 'PowerFreq_bkg_aas')
   
end