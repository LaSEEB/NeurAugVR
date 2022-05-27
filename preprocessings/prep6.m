function EEG = prep6(EEG,resamp,hp,lp,dirs,elims,ereject,no_interp_chans,no_discard_chans)
% E.g.:
% resamp = 250
% hp = 1
% lp = 40
% dirs = {'S  7','S  8'}  % Left and Right
% elims = [-5.5, 5.5] % (epoch limits [s], from arrow)
% ereject = true % true rejects epochs, false keeps them (and marks them)
% no_interp_chans = {'C3', 'C4'} (do not interpolate)
% no_discard_chans = {'C3', 'C4'} (do not discard to make data fully-rank)
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

%% Continuous clean
EEGtemp = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,'WindowCriterion',0.5,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 8] );
prep_report.('bursts') = EEG.xmax - EEGtemp.xmax; 
prep_report.('burstsP') = prep_report.('bursts')/EEG.xmax*100;
            
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
prep_report.('rej_comps') = numel(rej_vec);

%% Epoch
EEGdi = pop_epoch(EEG, dirs, elims, 'epochinfo', 'yes');
% Clean trials
typerej = 1;                % On EEG and not ICA components
elec_comp = 1:EEGdi.nbchan;
locthresh = 3;
globthresh = 3;
superpose = 1;              % Different than default
vistype = 0;
plotflag = 0;
[EEG, ~, ~, nrej1] = pop_jointprob(EEGdi, typerej, elec_comp, locthresh, globthresh, superpose, reject, vistype,[],plotflag);
[EEG, ~, ~, nrej2] = pop_rejkurt(EEG, typerej, elec_comp,locthresh, globthresh, superpose, ereject, vistype);

%% Report
urevents_before = [EEGdi.event(ismember({EEGdi.event.type},dirs)).urevent];
for di = 1:numel(dirs)
    urevent_di_before = [EEGdi.event(strcmp({EEGdi.event(:).type}, dirs{di})).urevent];
    urevent_di_after = [EEG.event(strcmp({EEG.event(:).type}, dirs{di})).urevent];
    kept = find(ismember(urevents_before,urevent_di_after));
    rej = find(ismember(urevents_before,urevent_di_before(~ismember(urevent_di_before,urevent_di_after))));
    prep_report.rej_trials.(regexprep(dirs{di}, ' ', '_')).kept = kept;
    prep_report.rej_trials.(regexprep(dirs{di}, ' ', '_')).rej = rej;
end

% for di = 1:numel(dirs)
%     prep_report.(strcat('trials', regexprep(dirs{di}, ' ', '_'))) = sum(strcmp({EEGdi.event(:).type}, dirs{di})) - sum(strcmp({EEG.event(:).type}, dirs{di}));
% end
% prep_report.('trials') = EEGdi.trials - EEG.trials;
% prep_report.('trialsP') = prep_report.('trials')/EEGdi.trials*100;
% prep_report.('total_trials') = EEGdi.trials;
% prep_report.('trials_rej_mask') = EEG.reject.rejjp | EEG.reject.rejkurt;

% fprintf(strcat('Prep 6 report\nChans removed: ',repmat('%s ',1,numel(prep_report.('chans'))),'\nBursts removed: %0.0f (%0.0f%%)\nComps removed: %d\nTrials removed: %0.0f / %d  (%0.0f%%)\n'),prep_report.('chans'){:},prep_report.('bursts'),prep_report.('burstsP'),prep_report.('comps'),prep_report.('trials'),prep_report.('total_trials'),prep_report.('trialsP'));
EEG.preproc = prep_report;

end


