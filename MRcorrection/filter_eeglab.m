function EEGFILT = filter_eeglab(EEG, low_edge, high_edge)
% Filter eeg with "low_edge" and "high_edge" using EEGLAB

    EEG_l = pop_eegfiltnew(EEG,low_edge,[]);
    EEGFILT = pop_eegfiltnew(EEG_l, [], high_edge);

end