function EEG = prep8(EEG,resamp,hp,lp)
% E.g.:
% resamp = 250
% hp = 1
% lp = 40
% Obs.: This preprocessing returns an epoched EEG!

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

%% Save keep channels
keep_chans = {'C3','C4'};
chns = find(ismember({EEG.chanlocs(:).labels}, keep_chans));
chans_data = [];
chans_locs = [];
for chi = 1:numel(keep_chans)
    chns(chi)
    chans_data = [chans_data; EEG.data(chns(chi),:)];
    chans_locs = [chans_locs; EEG.chanlocs(chns(chi))];
end

%% Remove bad channels
EEGallchans = EEG;
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',10,'ChannelCriterion',0.8,'LineNoiseCriterion',5,'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
% EEG = pop_select(EEG, 'nochannel',{'C3','EKG'}); % DEBUG
prep_report.('chans') = {EEGallchans.chanlocs(~ismember({EEGallchans.chanlocs(:).labels},{EEG.chanlocs(:).labels})).labels};

%% Add keep channels if necessary
for chi = 1:numel(keep_chans)
    if ~ismember(keep_chans{1},{EEG.chanlocs(:).labels})
        EEG.nbchan = EEG.nbchan+1;
        EEG.data(end+1,:) = chans_data(chi,:);
        EEG.chanlocs(1,EEG.nbchan)= chans_locs(chi);
    end
end
EEG = eeg_checkset(EEG);

%% Rank deficit
rank_deficit = EEGallchans.nbchan - EEG.nbchan;

%% Interpolate channels
EEG = pop_interp(EEG, EEGallchans.chanlocs, 'spherical');

%% Re-reference
EEG = fullRankAveRef(EEG);

%% Discard channels to make the data full ranked
if rank_deficit > 0
    keep_chans = {'C3','C4'};
    chns = find(ismember({EEG.chanlocs(:).labels}, keep_chans));
    channelSubset = loc_subsets(EEG.chanlocs, EEG.nbchan-rank_deficit,false,false,{chns});
    EEG = pop_select( EEG,'channel', channelSubset{1});
    EEG = pop_chanedit(EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
end

%% Continuous clean
EEGtemp = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,'WindowCriterion',0.5,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 8] );
prep_report.('bursts_bICA') = EEG.xmax - EEGtemp.xmax; 
prep_report.('burstsP_bICA') = prep_report.('bursts_bICA')/EEG.xmax*100;
            
%% ICA
EEGtemp = pop_runica(EEGtemp, 'icatype', 'runica','extended',1,'interrupt','on');

%% Weight transfer
EEG.icaweights = EEGtemp.icaweights;
EEG.icasphere = EEGtemp.icasphere;
EEG.icachansind = EEGtemp.icachansind;
EEG.icawinv = EEGtemp.icawinv;

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
prep_report.('comps') = numel(rej_vec);

%% Continuous clean: remove and interpolate bursts
EEGtemp = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,'WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
prep_report.('bursts_aICA_rej_mask') = sum(abs(EEG.data-EEGtemp.data),1) >= 1e-10;
prep_report.('bursts_aICA') = sum(prep_report.('bursts_aICA_rej_mask'))/EEG.srate;
prep_report.('burstsP_aICA') = prep_report.('bursts_aICA')/EEG.xmax*100;

%% Report
fprintf(strcat('Prep 8 report\nChans removed: ',repmat('%s ',1,numel(prep_report.('chans'))),'\nBursts removed before ICA: %0.0f (%0.0f%%)\nComps removed: %d\nBursts removed after ICA: %0.0f (%0.0f%%)\n'),prep_report.('chans'){:},prep_report.('bursts_bICA'),prep_report.('burstsP_bICA'),prep_report.('comps'),prep_report.('bursts_aICA'),prep_report.('burstsP_aICA'));
EEG.preproc = prep_report;

end