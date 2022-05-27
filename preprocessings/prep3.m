function EEG = prep3(EEG,resamp,hp,lp,no_interp_chans)
% E.g.:
% resamp = 250
% hp = 1
% lp = 40
% no_interp_chans = {'C3', 'C4'} (do not interpolate)

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

%% Save no-interp-chans channels
if isequal(no_interp_chans, 'all')
    no_interp_chans = {EEG.chanlocs(:).labels};
end

chns = find(ismember({EEG.chanlocs(:).labels}, no_interp_chans));
chans_data = [];
chans_locs = [];
for chi = 1:numel(no_interp_chans)
    chans_data = [chans_data; EEG.data(chns(chi),:)];
    chans_locs = [chans_locs; EEG.chanlocs(chns(chi))];
end

%% Remove bad channels
EEGallchans = EEG;
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',10,'ChannelCriterion',0.8,'LineNoiseCriterion',5,'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
% EEG = pop_select(EEG, 'nochannel',{'C3','EKG'}); % DEBUG

%% Add no-interp-channels if necessary
for chi = 1:numel(no_interp_chans)
    if ~ismember(no_interp_chans{1},{EEG.chanlocs(:).labels})
        EEG.nbchan = EEG.nbchan+1;
        EEG.data(end+1,:) = chans_data(chi,:);
        EEG.chanlocs(1,EEG.nbchan)= chans_locs(chi);
    end
end
EEG = eeg_checkset(EEG);
prep_report.('interp_chans') = {EEGallchans.chanlocs(~ismember({EEGallchans.chanlocs(:).labels},{EEG.chanlocs(:).labels})).labels};

%% Interpolate channels
EEG = pop_interp(EEG, EEGallchans.chanlocs, 'spherical');

%% Re-reference
EEG = fullRankAveRef(EEG);
% fprintf(strcat('Prep 3 report\nChans interpolated: ',repmat('%s ',1,numel(prep_report.('interp_chans'))),'\n'),prep_report.('interp_chans'){:});
EEG.preproc = prep_report;

end


