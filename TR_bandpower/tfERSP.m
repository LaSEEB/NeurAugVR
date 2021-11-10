function [tr_vecs, ersp_vecs] = tfERSP(EEG, chan, vol_event, bands, trial_ntr, trial_tr_id)
% E.g.:
% chan = 'C3';
% bands = [1, 4; 4, 8; 8, 13; 13, 30];
% vol_event = 'Vol Start';
% trial_ntr = 10;  % Number of TRs to calculate a ERSP window with
% trial_tr_id = 8; % In that number of TRs, TR that should be attributed to the resulting ERSP value

% Channel
chn = find(ismember({EEG.chanlocs(:).labels}, chan));

% Volume latencies
vol_lats = [EEG.event(strcmp({EEG.event(:).type},vol_event)).latency];
nvols = numel(vol_lats);

tr = 1.26;
fs = EEG.srate;
tr_len = tr * fs;
vol_lats = [vol_lats, vol_lats(end) + tr_len];  % Add an extra latency at the end, to calculate the power of the last volume

base_ntr = floor(trial_ntr/2); % the rest of trial-tr's are considered task
base_len = base_ntr * tr_len;
% tr_vec = zeros(1,nvols-trial_ntr+1);

% Power
ersp_vecs = zeros(size(bands, 1),nvols-trial_ntr+1);

% TF settings
baseline = nan; % Indiferent actually
scale = 'abs';
cycles = 0;
ylims = [0, 30];
        
for v = 1:numel(vol_lats) - trial_ntr
    %% Get trial
    start = vol_lats(v);
    finish = vol_lats(v+trial_ntr) - 1;
    if v == numel(vol_lats) - trial_ntr
        lala=0;
    end
    tubmat = double(EEG.data(chn, start:finish)); % I use the word 'tub' for a 1x1xN array, and 'tubmat' for a 1xMxN array
    timelims = [EEG.times(start) EEG.times(finish)] - EEG.times(start+base_len);
    
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
        ersp_vecs(b,v) = bandpow_task/bandpow_base*100;
    end
%     tr_vec(v) = v - 1 + trial_tr_id;
    % tr_vec = [tr_vec, v - 1 + trial_tr_id];
    
end

% Add missing start- and finish- TRs with 0 ERDs
ersp_vecs = [zeros(size(bands,1),trial_tr_id-1), ersp_vecs, zeros(size(bands,1),trial_ntr-trial_tr_id)];
% Add times, with 0 in first TR
tr_vecs(1,:) = 1:nvols;
tr_vecs(2,:) = (EEG.times(vol_lats(tr_vecs(1,:))) - EEG.times(vol_lats(1)))/1000;

% Normalize ERSP vecs
for b = 1:size(bands,1)
    ersp_vecs(b,:) = -1 + 2.*(ersp_vecs(b,:) - min(ersp_vecs(b,:)))./(max(ersp_vecs(b,:)) - min(ersp_vecs(b,:)));
end

end