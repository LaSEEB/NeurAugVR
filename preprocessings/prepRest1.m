function EEG = prepRest1(EEG,resamp,hp,lp,no_interp_chans,no_discard_chans)
% E.g. 1:
% resamp = 250
% hp = 1
% lp = 40
% no_interp_chans = {'C3', 'C4'} (do not interpolate)
% no_discard_chans = {'C3', 'C4'} (do not discard to make data fully-rank)

% E.g. 2:
% no_interp_chans = {}
% no_discard_chans = 'all'

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

%% Rank deficit
rank_deficit = EEGallchans.nbchan - EEG.nbchan;

%% Interpolate channels
EEG = pop_interp(EEG, EEGallchans.chanlocs, 'spherical');

%% Re-reference
EEG = fullRankAveRef(EEG);

%% Discard channels to make the data full ranked
EEGallchans = EEG;
if rank_deficit > 0
    if isequal(no_discard_chans, 'all')
        no_discard_chans = {EEG.chanlocs(:).labels};
    end
    chns = find(ismember({EEG.chanlocs(:).labels}, no_discard_chans));
    channelSubset = loc_subsets(EEG.chanlocs, EEG.nbchan-rank_deficit,false,false,{chns});
    EEG = pop_select( EEG,'channel', channelSubset{1});
    EEG = pop_chanedit(EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
end
prep_report.('rej_chans') = {EEGallchans.chanlocs(~ismember({EEGallchans.chanlocs(:).labels},{EEG.chanlocs(:).labels})).labels};

%% ICA
EEG = pop_runica(EEG, 'icatype', 'runica','extended',1,'interrupt','on');

%% Prun
EEG = iclabel(EEG);
iclabel_mat = EEG.etc.ic_classification.ICLabel.classifications; % n*7 matrix where n=number of IC's and 7=number o classes
thres = 0.85; % 90% o.o
rej_vec = [];
for n = 1:size(iclabel_mat,1) % For each IC, determine if it is to reject
    [val, idx] = max(iclabel_mat(n, :));
    if val >= thres && (idx == 2 || idx == 3)
        rej_vec = [rej_vec, n];
    end
end
EEG = pop_subcomp(EEG, rej_vec, 0);
prep_report.('rej_comps') = numel(rej_vec);

%% Continuous clean: remove and interpolate bursts
EEGtemp = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,'WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
bursts_aICA_rej_mask = sum(abs(EEG.data-EEGtemp.data),1) >= 1e-10;
bursts_aICA = sum(bursts_aICA_rej_mask)/EEG.srate;
prep_report.('rej_segments') = bursts_aICA/EEG.xmax*100;
prep_report.('rej_segments_mask') = bursts_aICA_rej_mask;

EEG = EEGtemp;
% - If you want to interpolate bad segments, use: 'BurstRejection','off' [default]
% - If you want to remove bad segments, use: 'BurstRejection','on'
% - If you don't want to remove/interpolate bad segments, but just identify
% them, comment the previous line "EEG = EEGtemp", and after preprocessing,
% check EEG.preproc.rej_segments_mask to know which time instances were
% interpolated (1) and which were not (0)

%% Report
EEG.preproc = prep_report;

end
