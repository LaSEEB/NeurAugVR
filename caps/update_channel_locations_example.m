%% Add
addpath('C:/Users/guta_/Documents/MATLAB/eeglab2021.0')  % Change according to the computer
varsbefore = who; eeglab; close; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})

%% Load
chanlocs_path = 'chanlocsMR32.mat';
data_path = 'EEG.mat';  % Just as an example: this file does not exist on github
load(chanlocs_path,'chanlocs');
load(data_path,'EEG');

%% Update
if isempty(EEG.chanlocs)  % Insert chanloc struct
    EEG.chanlocs = chanlocs;
else  % Insert chanloc struct matching the label order of already existing EEG.chanloc
    chanlocs_temp = [];
    for chi = 1:numel(EEG.chanlocs)
        chanlocs_temp = [chanlocs_temp, chanlocs(ismember({chanlocs(:).labels}, EEG.chanlocs(chi).labels))];
    end
    EEG.chanlocs = chanlocs_temp;
end

%% Check
EEG = eeg_checkset(EEG);
