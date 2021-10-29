function EEG = prep5(EEG,resamp,hp,lp,dirs,elims)
% E.g.:
% resamp = 250
% hp = 1
% lp = 40
% dirs = {'S  7','S  8'}  % Left and Right
% elims = [-5.5, 5.5] % (epoch limits [s], from arrow)
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

%% Report
prep_report.('comps') = numel(rej_vec);

%% Epoch
EEGdi = pop_epoch(EEG, dirs, elims, 'epochinfo', 'yes');

%% Clean trials
typerej = 1;                % On EEG and not ICA components
elec_comp = 1:EEGdi.nbchan;
locthresh = 3;
globthresh = 3;
superpose = 1;              % Different than default
reject = 1;
vistype = 0;
plotflag = 0;
[EEG, ~, ~, nrej1] = pop_jointprob(EEGdi, typerej, elec_comp, locthresh, globthresh, superpose, reject, vistype,[],plotflag);
[EEG, ~, ~, nrej2] = pop_rejkurt(EEG, typerej, elec_comp,locthresh, globthresh, superpose, reject, vistype);

%% Report
for di = 1:numel(dirs)
    prep_report.(strcat('trials', regexprep(dirs{di}, ' ', '_'))) = sum(strcmp({EEGdi.event(:).type}, dirs{di})) - sum(strcmp({EEG.event(:).type}, dirs{di}));
end
prep_report.('trials') = EEGdi.trials - EEG.trials;
prep_report.('trialsP') = prep_report.('trials')/EEGdi.trials*100;
prep_report.('total_trials') = EEGdi.trials;

fprintf(strcat('Prep 5 report\nChans removed: ',repmat('%s ',1,numel(prep_report.('chans'))),'\nComps removed: %d\nTrials removed: %0.0f / %d  (%0.0f%%)\n'),prep_report.('chans'){:},prep_report.('comps'),prep_report.('trials'),prep_report.('total_trials'),prep_report.('trialsP'));
EEG.preproc = prep_report;

end


