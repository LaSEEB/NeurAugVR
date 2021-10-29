function [tr_vecs, pow_vecs] = tfiltpower(EEG, chan, vol_event, bands)
% Channel
chn = find(ismember({EEG.chanlocs(:).labels}, chan));

% Volume latencies
vol_lats = [EEG.event(strcmp({EEG.event(:).type},vol_event)).latency];
nvols = numel(vol_lats);
avg_lat_dif = sum(diff(vol_lats))/(nvols-1);
vol_lats = [vol_lats, vol_lats(end) + avg_lat_dif];  % Add an extra latency at the end, to calculate the power of the last volume

% Power
pow_vecs = zeros(size(bands, 1),nvols);

% TR
tr_vecs = zeros(2,nvols);  % Indexes and times
tr_vecs(1,:) = 1:nvols;
tr_vecs(2,:) = (EEG.times(vol_lats(tr_vecs(1,:))) - EEG.times(vol_lats(1)))/1000;

% Calculate
for b = 1:size(bands, 1)
    EEGhp = pop_eegfiltnew(EEG, 'locutoff',bands(b,1), 'plotfreqz',0);
    EEGbp = pop_eegfiltnew(EEGhp, 'hicutoff',bands(b,2), 'plotfreqz',0);
    
    for v = 1:numel(vol_lats)-1
        pow_vecs(b,v) = mean(EEGbp.data(chn,vol_lats(v):vol_lats(v+1)).^2);
    end
    pow_vecs(b,:) = -1 + 2.*(pow_vecs(b,:) - min(pow_vecs(b,:)))./(max(pow_vecs(b,:)) - min(pow_vecs(b,:)));
    
end
end