%% Add paths
clear, clc
addpath('C:/Users/guta_/Documents/MATLAB/eeglab2021.0')
varsbefore = who; eeglab; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})

%% Load
EEG = pop_loadbv('..\data\', 'sub-pilot03_ses-outside_task-neurowMIMO.vhdr', [], []);

% If necessary, take out auxiliary channel
% aux_start = 33;
% EEG.dataaux = EEG.data(aux_start:end,:);
% EEG.data = EEG.data(1:aux_start-1,:);
% EEG.chanlocsaux = EEG.chanlocs(aux_start:end);
% EEG.chanlocs = EEG.chanlocs(1:aux_start-1);
% EEG.nbchan = size(EEG.data, 1);

%% Preprocess (prep. 1)
bandlims = [1, 40];
resamp = 250;
EEG = pop_select(EEG, 'nochannel',{'ECG','EKG'});
EEG = pop_resample(EEG, resamp);
EEG = pop_eegfiltnew(EEG, 'locutoff',bandlims(1), 'plotfreqz',0);
EEG = pop_eegfiltnew(EEG, 'hicutoff',bandlims(2), 'plotfreqz',0);
EEG = pop_reref(EEG, []);

%% Calculate/Plot topo-ERSP%
dirs = {'S  7', 'S  8'};
figure
th = tiledlayout(1,numel(dirs));
plot_options = struct('label',{'ref','gnd','c3','c4'},'show_label',{1,1,1,1},'marker',{'x','x','.','.'},'marker_col',{[1,1,1], [1,1,1], [0,0,0], [0,0,0]});

for di = 1:numel(dirs)
    %% Initialize topoplot parameters
    erd_array = zeros(1,numel(EEG.chanlocs));
    loc_array = [];
    
    %% Epoch
    dire = dirs(di);  % Left
    epochlims = [-5.5, 5.5];
    EEGdi = pop_epoch(EEG, dire, epochlims, 'epochinfo', 'yes');
    
    for chn = 1:numel(EEG.chanlocs)
        %% Select channel
        tubmat = double(EEGdi.data(chn, :, :));  % I use the word 'tub' for a 1x1xN array, and 'tubmat' for a 1xMxN array
        
        %% Time-Frequency decomposition of each trial
        baseline = nan; % Indiferent actually
        scale = 'abs';
        cycles = 0;
        [~, ~, ~, time_vec, freq_vec, ~, ~, amp_vol] = newtimef_trueamp(tubmat, size(tubmat,2), [EEGdi.times(1) EEGdi.times(end)], EEGdi.srate, cycles,'plotitc','off', 'baseline', baseline, 'scale', scale, 'plotersp', 'off');
        time_vec = time_vec/1000;
        amp_vol = double(amp_vol);
        
        %% Calculate power
        pow_vol = amp_vol.^2;
        
        %% Average trials
        pow_mat = mean(pow_vol,3);
        
        %% Calculate ERD%
        baselims = [-4, 0];
        erd_mat = (pow_mat./mean(pow_mat(:,time_vec > baselims(1) & time_vec < baselims(2)),2) - 1)*100;
        
        %% Average alpha's frequency band
        freqlims = [8, 12];
        erd_vec = mean(erd_mat(freq_vec > freqlims(1) & freq_vec < freqlims(2),:),1);

        %% Average task's time interval
        tasklims = [0.5, 4];
        erd_val = mean(erd_vec(time_vec > tasklims(1) & time_vec < tasklims(2)));
        
        %% Get channel location
        loc_struct = EEGdi.chanlocs(chn);
        
        %% Set default plot options
        loc_struct.show_label = 0;
        loc_struct.marker = '.';
        loc_struct.marker_col = [0,0,0];
        
        %% Set specified plot options
        for opt = 1:numel(plot_options)
            if strcmpi(strtrim(loc_struct.labels),plot_options(opt).label)
                loc_struct.show_label = plot_options(opt).show_label;
                loc_struct.marker = plot_options(opt).marker;
                loc_struct.marker_col = plot_options(opt).marker_col;
            end
        end
        
        %% Store
        erd_array(chn) = erd_val;
        loc_array = [loc_array; loc_struct];
        
    end
    %% Plot
    nexttile
    topoplot_custom_chans(erd_array, loc_array, 'electrodes', 'labelpoint');
    % Optionally add colored title
    title(sprintf('\\bf{%s}\\color{%s}{\\bf{%s}}', dire{1}(1:end-1), 'blue', dire{1}(end)), 'interpreter', 'tex')
end
% Optionally save space
% th.TileSpacing = 'none';
