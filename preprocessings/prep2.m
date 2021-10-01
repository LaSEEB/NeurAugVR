function EEG = prep2(EEG,resamp,hp,lp)
% E.g.:
% resamp = 250
% hp = 1
% lp = 40

%% Remove ECG
EEG = pop_select(EEG, 'nochannel',{'ECG','EKG'});

%% Resample
if ~isempty(resamp)
    if EEG.srate ~= resamp
        EEG = pop_resample(EEG, resamp);
    end
end

%% Filter
EEG = pop_eegfiltnew(EEG, 'locutoff',hp, 'plotfreqz',0);
EEG = pop_eegfiltnew(EEG, 'hicutoff',lp, 'plotfreqz',0);

%% Re-reference
EEG = fullRankAveRef(EEG);
EEG.preproc = [];

end


