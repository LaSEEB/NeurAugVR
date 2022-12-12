function wins = segment_signal(eeg,ecg_markers_win,win_max_len)

wins = nan(numel(ecg_markers_win)-1,win_max_len);

for i = 1:numel(ecg_markers_win)-1
  start = ecg_markers_win(i);
  finish = ecg_markers_win(i+1);
  segment = eeg(start:finish);
  finish_trimmed = min(numel(segment),win_max_len);
  wins(i,1:finish_trimmed) =   segment(1:finish_trimmed);
end
end
