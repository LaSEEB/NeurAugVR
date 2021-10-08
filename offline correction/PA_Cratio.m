function opt_param = PA_Cratio(filelocation, filename, projic_k, obs_k, obs_npc, aas_k, aas_win, plt) 

% Ballistocardiogram artifact correction taking into account physiological
% signal preservation in simultaneous EEG-fMRI
% Abreu et al., 2016
% NeuroImage 135 (2016) 4563
%%
% PowerFreq(1, :) - power - (N channels x art_harm cardiac harmonics)
% PowerFreq(2, :) - params
% PowerFreq(3, :) - ICs


    %% PROJIC  

    load([filelocation, '/PowerFreq_bkg_',filename])
    load([filelocation, '/PowerFreq_', filename])

    nr_k = length(projic_k);
    k_vec =  2:length(PowerFreq);

    r_bkg = zeros(1, length(PowerFreq)-1);
    r_art = zeros(1, length(PowerFreq)-1);

    % uncorrected corresponds to the first element (k = [])
    s_unc_art = sum(PowerFreq{1, 1}(:));
    s_unc_bkg = sum(PowerFreq_bkg{1, 1}(:));

    % for each k (number of clusters) 
    for k =  k_vec

        s_corr_art = sum(PowerFreq{1, k}(:));
        s_corr_bkg = sum(PowerFreq_bkg{1, k}(:));
        r_art(1, k-1) = (s_unc_art - s_corr_art)/s_unc_art; % artifact ratio
        r_bkg(1, k-1) = (s_unc_bkg - s_corr_bkg)/s_unc_bkg; % physiological ratio
    end

    % varying w_bkg (the weight for the physiological signal):
    % 0 - the BCG artifact removal is the priority
    % 1 - the physiological background signal preservation is the priority 
    w_vec = 0:0.1:1;
    c_vec = zeros(length(w_vec), nr_k);
    for w = 1:length(w_vec)
        w_bkg = w_vec(w);
        c = w_bkg*(1-r_bkg) + (1-w_bkg)*r_art; % combined ratio C(w_bkg)
        c_vec(w, :) = c; % C values for each weight w_bkg
    end

    % optimization:
    [val, k_ind] = max(c_vec, [], 2);

    % maximum cratio for each weight:
    cratio_projic = val;

    % optimal k's for each weight:
    opt_k_vec = k_vec(k_ind);

    if plt == 1
        st = cellfun( @num2str, num2cell(projic_k), 'UniformOutput',false);
        C = cell(1, length(k_vec));
        C(:) = {'k = '};
        leg = strcat(C, st);

        figure('Position', get(0, 'Screensize'));
        bar(0:0.1:1, c_vec)
        legend(leg)
        ylabel('Combined ratio')
        ylim([0 1])
        xlabel('Background weight')
        title(['Parameter optimization - PROJIC - ', filename])
        print([path, '/Figures/', filename, '_ratioC_projic_mod'], '-dpng')
    end

    %% 'PROJIC-OBS'

    load([filelocation, '/PowerFreq_bkg_obs_', filename])
    load([filelocation, '/PowerFreq_obs_', filename])

    nr_k = length(obs_k);
    nr_pc = length(obs_npc);

    %%%%%%%
    r_bkg = zeros(nr_k, nr_pc);
    r_art = zeros(nr_k, nr_pc);

    % uncorrected corresponds to the first element (k = [], no OBS)
    s_unc_art = sum(PowerFreq_obs{1, 1, 1}(:));
    s_unc_bkg = sum(PowerFreq_bkg_obs{1, 1, 1}(:));


    for k = 2:nr_k+1
        for n = 1:nr_pc
            s_corr_art = sum(PowerFreq_obs{1, n, k}(:));
            s_corr_bkg = sum(PowerFreq_bkg_obs{1, n, k}(:));
            r_art(k-1, n) = (s_unc_art - s_corr_art)/s_unc_art;
            r_bkg(k-1, n) = (s_unc_bkg - s_corr_bkg)/s_unc_bkg;
        end
    end


    w_vec = 0:0.1:1;
    c_vec = zeros(length(w_vec), nr_k*nr_pc);
    cratio_projic_obs = zeros(length(w_vec), 1);
    opt_param_obs = zeros(length(w_vec), 2);
    for w = 1:length(w_vec)
        w_bkg = w_vec(w);
        c = (w_bkg*(1-r_bkg) + (1-w_bkg)*r_art)';
        c_vec(w, :) = c(:);
        [val, ind] = max(c(:));
        [row, col] = ind2sub(size(c),ind);
        cratio_projic_obs(w) = val;
        opt_param_obs(w, :) = [row, col]; % optimal nr_pc and k for each weight
    end

    if plt == 1
        %%% plot legend
        st_k = cellfun(@num2str, num2cell(reshape(repmat(obs_k, nr_pc, 1), 1, nr_pc*nr_k)), 'UniformOutput',false);
        C_k = cell(1, nr_k*nr_pc);
        C_k(:) = {'k = '};

        st_pc = cellfun( @num2str, num2cell(reshape(repmat(obs_npc, 1, nr_k), 1, nr_pc*nr_k)), 'UniformOutput',false);
        C_pc = cell(1, nr_k*nr_pc);
        C_pc(:) = {', pc = '};
        leg_obs = strcat(C_k, st_k, C_pc, st_pc);
        %%%

        figure('Position', get(0, 'Screensize'));
        bar(0:0.1:1, c_vec)
        legend(leg_obs);
        ylabel('Combined ratio')
        ylim([0 1])
        xlabel('Background weight')
        title(['Parameter optimization - PROJIC-OBS - ', filename])
        print([path, '/Figures/', filename, '_ratioC_projic-obs'], '-dpng')

    end
    %% 'PROJIC-AAS'
    load([filelocation, '/PowerFreq_bkg_aas_', filename])
    load([filelocation, '/PowerFreq_aas_', filename])


    nr_k = length(aas_k);
    nr_w = length(aas_win);
    r_bkg = zeros(nr_k, nr_w);
    r_art = zeros(nr_k, nr_w);

    s_unc_art = sum(PowerFreq_aas{1, 1, 1}(:));
    s_unc_bkg = sum(PowerFreq_bkg_aas{1, 1, 1}(:));


    for k =  2:length(PowerFreq_aas)
        for n = 1:nr_w
            s_corr_art = sum(PowerFreq_aas{1, n, k}(:));
            s_corr_bkg = sum(PowerFreq_bkg_aas{1, n, k}(:));
            r_art(k-1, n) = (s_unc_art - s_corr_art)/s_unc_art;
            r_bkg(k-1, n) = (s_unc_bkg - s_corr_bkg)/s_unc_bkg;
        end
    end


    w_vec = 0:0.1:1;
    c_vec = zeros(length(w_vec), nr_k*nr_w);
    opt_param_aas = zeros(length(w_vec), 2);
    cratio_projic_aas = zeros(length(w_vec), 1);
    for w = 1:length(w_vec)
        w_bkg = w_vec(w);
        c = (w_bkg*(1-r_bkg) + (1-w_bkg)*r_art)';
        c_vec(w, :) = c(:);
        [val, ind] = max(c(:));
        cratio_projic_aas(w) = val;
        [row, col] = ind2sub(size(c),ind);
        opt_param_aas(w, :) = [row, col]; % optimal nr_w and k for each weight
    end


    if plt == 1
        %%% plot legend
        st_k = cellfun(@num2str, num2cell(reshape(repmat(aas_k, nr_w, 1), 1, nr_w*nr_k)), 'UniformOutput',false);
        C_k = cell(1, nr_w*nr_k);
        C_k(:) = {'k = '};

        st_w = cellfun( @num2str, num2cell(reshape(repmat(aas_win, 1, nr_k), 1, nr_w*nr_k)), 'UniformOutput',false);
        C_w = cell(1, nr_w*nr_k);
        C_w(:) = {', win = '};
        leg_aas = strcat(C_k, st_k, C_w, st_w);
        %%%

        figure('Position', get(0, 'Screensize'));
        bar(0:0.1:1, c_vec)
        legend(leg_aas);
        ylabel('Combined ratio')
        ylim([0 1])
        xlim([-0.1 1.1])
        xlabel('Background weight')
        title(['Parameter optimization - PROJIC-AAS - ', filename])
        print([path, 'Figures/' , filename, '_ratioC_projic-aas'])
    end 
    %%
    k = opt_k_vec';
    c = cratio_projic;

    AAS_k = aas_k(opt_param_aas(:, 2))';
    AAS_win = aas_win(opt_param_aas(:, 1))';
    AAS_c = cratio_projic_aas;

    OBS_k = obs_k(opt_param_obs(:, 2))';
    OBS_pc = obs_npc(opt_param_obs(:, 1))';
    OBS_c = cratio_projic_obs;

    weight = w_vec';
    opt_param = table(weight, k, c, OBS_k, OBS_pc, OBS_c, AAS_k, AAS_win, AAS_c);    


        
