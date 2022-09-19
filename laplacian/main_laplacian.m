%% Add EEGLAB path
addpath('C:/Users/guta_/Documents/MATLAB/eeglab2021.0');
varsbefore = who; eeglab; close; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})

%% Load data
load('EEG.mat','EEG')

%% Extract XYZ coordinates from EEG structure
X = [EEG.chanlocs.X];
Y = [EEG.chanlocs.Y];
Z = [EEG.chanlocs.Z];

%% Calculate Laplacian
[surf_lap,~,~] = laplacian_perrinX(EEG.data(:,:,1),X,Y,Z,[],1e-5);

%% Plot
figure
topoplot(surf_lap,EEG.chanlocs,'plotrad',.53,'maplimits',[-40 40],'electrodes','off');
