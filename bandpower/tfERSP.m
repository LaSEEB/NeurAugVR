function [x_vecs, ersp_vecs] = tfERSP(EEG, chan, onsets, bands, trial, time_id)
%% Data settings
chn = find(ismember({EEG.chanlocs(:).labels}, chan));
data = EEG.data;
fs = EEG.srate;

%% Trial settings
base_len = abs(trial(1)) * fs;
task_len = trial(2) * fs;
id = round(time_id * fs) + 1;

%% Power
x_vecs = zeros(2,numel(onsets));
x_vecs(1,:) = 1:numel(onsets);          % Onset ids
x_vecs(2,:) = EEG.times(onsets)/1000;   % Onset times
x_vecs(3,:) = onsets;                   % Onset latencies
ersp_vecs = zeros(size(bands, 1),numel(onsets));

%% TF settings
baseline = nan; % Indiferent actually
scale = 'abs';
cycles = 0;
ylims = [0, 30];

%% Slide
for n = 1:numel(onsets)
    
    start = onsets(n) - base_len;
    finish = onsets(n) + task_len - 1;
    
    if start > 0 && finish < size(data, 2) + 1 
        
        %% Calculate shift, to save the ERSP value in onset + shift
        onsets_in_win = onsets(onsets >= start & onsets <= finish) - onsets(n);
        shift = find(onsets_in_win < id, 1, 'last') - numel(onsets(onsets >= start & onsets <= onsets(n)));
        
        %% Get data and times
        tubmat = double(EEG.data(chn, start:finish)); % I use the word 'tub' for a 1x1xN array, and 'tubmat' for a 1xMxN array
        timelims = [EEG.times(start) EEG.times(finish)] - EEG.times(onsets(n));

        %% Time-Frequency decomposition of each trial
        [~, ~, ~, time_vec, freq_vec, ~, ~, amp_mat] = newtimef_trueamp(tubmat, size(tubmat,2), timelims, EEG.srate, cycles, 'freqs', ylims,'plotitc','off', 'baseline', baseline, 'scale', scale, 'plotersp', 'off');
        time_vec = time_vec/1000;
        amp_mat = double(amp_mat);
        
        %% Calculate power
        pow_mat = amp_mat.^2;
        
        %% Calculate ERD%
        for b = 1:size(bands,1)
            baselims = [-4, 0];
            tasklims = [0.5, 4];
            bandpow_vec = mean(pow_mat(freq_vec > bands(b,1) & freq_vec < bands(b,2),:),1);
            bandpow_base = mean(bandpow_vec(time_vec > baselims(1) & time_vec < baselims(2)),2);
            bandpow_task = mean(bandpow_vec(time_vec > tasklims(1) & time_vec < tasklims(2)),2);
            ersp_vecs(b,n + shift) = bandpow_task/bandpow_base*100;
        end
    end
end

%% Normalize ERSP vecs
for b = 1:size(bands,1)
    ersp_vecs(b,:) = -1 + 2.*(ersp_vecs(b,:) - min(ersp_vecs(b,:)))./(max(ersp_vecs(b,:)) - min(ersp_vecs(b,:)));
end
end