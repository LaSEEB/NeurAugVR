addpath('/home/iesteves/MR_Correction_Pipeline/eeglab2019_0') 
eeglab

%%  PA Correction using optimal parameters

%cd '/home/iesteves/MR_Correction_Pipeline/Preprocessing/PAfiltered_PROJIC_OBS/'
%files = dir('*KDEF*_1*.set');
cd '/home/iesteves/MR_Correction_Pipeline/Preprocessing/Bireref/'
files = dir('*.mat');

filenames = {files.name};
filefolder = {files.folder};

path = '/home/iesteves/MR_Correction_Pipeline/Preprocessing/';

locfilepath = '/home/iesteves/MR_Correction_Pipeline/Pipeline';
locfilename = strcat(locfilepath, '/chan_loc.ced'); %'Chan_loc_x20t.ced';

chan = 6; % change channel number for plots
ecgchan = 32;

%% EEG filter 0.3-70Hz + bireref + ICA weights
nchan  = 31;
for n = 1:length(filenames)
    
    strfile = strsplit(filenames{n}, {'_', '.'});
    subject = strfile{1};
    session = strfile{2};
    task = strfile{3};
    
    
    EEG = pop_loadset('filename',filenames{n}, 'filepath', filefolder{n});
    
    
    EEG = filter_eeglab(EEG, 0.3, 70);
    EEG.bpfilter = [0.3  70];

    ecgchan = 32;
    el = EEG.data([1:ecgchan-1], :);
    [ bi_mean, ~] = myBiweight(el');

    EEG.data([1:ecgchan-1], :) = el - repmat(bi_mean, nchan, 1); % re-referencing
    EEG.bimean = bi_mean;
    
    EEG = pop_runica(EEG, 'icatype', 'runica', 'chanind', [1:31], 'extended', 1);
    save([path, '/Bireref/',subject, '_', session, '_',  task, '_filtered_bireref', '.mat'], 'EEG')

end
%%

%% IClabel IC rejection
nchan  = 31;
channels = 1:31;
for k = 1:length(filenames)
    
    strfile = strsplit(filenames{k}, {'_', '.'});
    subject = strfile{1};
    session = strfile{2};
    task = strfile{3};
    
    load([path, '/Bireref/', subject, '_', session, '_',  task, '_filtered_bireref', '.mat'])
    
    EEG.data = double(EEG.data);
    % icalabel
    EEG = iclabel(EEG);

    % classifications
    EEG.etc.ic_classification.ICLabel.classifications % n*7 matrix where n=number of IC's and 7=number o classes

    % show classes
    EEG.etc.ic_classification.ICLabel.classes
    % ans = 'Brain'    'Muscle'    'Eye'    'Heart'    'Line Noise'    'Channel Noise'    'Other'

    %IC threshold value in %
    thres = 0.9; % 50%
    
    % Find IC's to reject
    for n = 1:length(EEG.etc.ic_classification.ICLabel.classifications)

        % find the index of the maximum classification score for each IC
        [val, idx] = max(EEG.etc.ic_classification.ICLabel.classifications(n, :));
        % if val > 50% and idx not 1 (Brain), then add to list.
        % todo:check also brain percentage
        if (val > thres) && (all(idx ~= [1, 7]))
            % add to a rejection list
            rej_ic(:,n) = n; % add the IC number
        else
            rej_ic(:,n) = NaN;
        end

    end

    rej_ic = rej_ic(~isnan(rej_ic)); %clear the NaNs
    
    % reject components
    EEG = pop_subcomp(EEG, rej_ic, 0);
    EEG.rej_ic = rej_ic;
    save([path, '/ICremoved/',subject, '_', session, '_',  task, '_ICremoved', '.mat'], 'EEG')
    
    signal = EEG.data(channels, :);
    [newsig_capmean, ~, ~] = outlier_rejection(signal, 'capmean');

    EEGCAP = EEG; EEGCAP.data(channels, :) = newsig_capmean;
    EEGCAP.outlier.method = 'cap';
    EEGCAP.outlier.win = '1';
    EEG = EEGCAP;
    save([path, '/Preproc/',subject, '_', session, '_',  task, '_preproc', '.mat'], 'EEG')

end