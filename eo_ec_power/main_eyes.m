%% Add 
clear, clc, close all
restoredefaultpath; clear RESTOREDEFAULTPATH_EXECUTED;  % Remove any added paths
addpath('Power spectrum')
addpath('C:/Users/guta_/Documents/MATLAB/eeglab2021.0')
varsbefore = who; eeglab; close; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})

%% Load
fol = 'C:\Users\guta_\Desktop\Data Analysis\Data\NeurAugVR\Formatted\';
fil = 'sub-05_ses-01_run-01';
load(strcat(fol,fil),'EEG');

%% Resample
fs = 250;
if EEG.srate ~= fs
    EEG = pop_resample(EEG, fs);
    for m = 1:size(EEG.event,2)
        EEG.event(m).latency = round(EEG.event(m).latency);
    end
end

%% Prepare EO/EC markers for epoching for Inês
marker = '3385';
epoch_dur = 1*60;
c = 1;
for i = 1:numel(EEG.event)
    if strcmp(EEG.event(i).type,marker)
        EEG.event(i).duration = epoch_dur;
        switch c
            case 1, EEG.event(i).label = 'EO';
            case 2, EEG.event(i).label = 'EC';
            case 3, EEG.event(i).label = 'EO';
            case 4, EEG.event(i).label = 'EC';
        end 
        c = c + 1;
    end
end

% % Add/update EO markers
% start = 1;
% for i = 1:numel(ec_ids)
%     
%     EEG.event(ec_ids(i)).duration = epoch_dur;
%     EEG.event(ec_ids(i)).label = 'EC';
%     
%     eo_ids = find(ismember({EEG.event(start:ec_ids(i)).type},'Open eyes'));
%     
%     switch numel(eo_ids)
%         case 0
%             lat = EEG.event(ec_ids(i)).latency - epoch_dur*EEG.srate;
%             EEG = add_event(EEG,'Open eyes',lat);
%             EEG.event(ismember({EEG.event(start:ec_ids(i)).type},'Open eyes')).duration = epoch_dur;
%             EEG.event(ismember({EEG.event(start:ec_ids(i)).type},'Open eyes')).label = 'EO';
%         case 1
%             EEG.event(eo_ids(1)).duration = epoch_dur;
%             EEG.event(eo_ids(1)).label = 'EO';
%         otherwise
%             disp('Too many Open eyes markers found!')
%             break
%     end
%     start = ec_ids(i) + 1;
% end
% 
% % Remove extra EO markers at the end
% EEG.event(ismember({EEG.event(start:end).type},'Open eyes')) = [];
% for i = 1:size(EEG.event,1)
%     EEG.event(i).urevent = i;
% end

%% Epoch
EEGEO = epoch_and_merge(EEG, 'EO');
EEGEC = epoch_and_merge(EEG, 'EC');

%% Calculate power
dur = 90;
chns = 20;  % 9-O1, 10-O2, 20-Oz;
% chns = 1:numel(EEG.chanlocs)-1;
method = 'jsousa';

fEO_mat = [];
fEC_mat = [];
pxxEO_mat = [];
pxxEC_mat = [];
f_alpha_mat = [];
p_alpha_mat = [];
alpha_ratio_mat = [];

start = 20*EEG.srate;
finish = start + dur*EEG.srate - 1;

for chn = chns
    [fEO, fEC, pxxEO, pxxEC, f_alpha, p_alpha, band_alpha] = iab_features(EEGEO.data(chn,start:finish), EEGEC.data(chn,start:finish), EEG.srate, method);
    alpha_ratio = (band_alpha(1) - band_alpha(2)) ./ band_alpha(2) * 100;
    
    fEO_mat = [fEO_mat, fEO];
    fEC_mat = [fEC_mat, fEC];
    pxxEO_mat = [pxxEO_mat, pxxEO];
    pxxEC_mat = [pxxEC_mat, pxxEC];
    f_alpha_mat = [f_alpha_mat, f_alpha];
    p_alpha_mat = [p_alpha_mat, p_alpha];
    alpha_ratio_mat = [alpha_ratio_mat, alpha_ratio];
end

%% Plot
chn_id = 1; % Id of the chan to plot from the previously selected chans (if all, chn_id = 20 is Oz)
figure
f_alpha = f_alpha_mat(:,chn_id);
p_alpha = p_alpha_mat(:,chn_id);
fEO = fEO_mat(:,chn_id);
pxxEO = pxxEO_mat(:,chn_id);
fEC = fEC_mat(:,chn_id);
pxxEC = pxxEC_mat(:,chn_id);

% nexttile
ltf = round(f_alpha(1), 2);
iaf = round(f_alpha(2), 2);
htf = round(f_alpha(3), 2);

ymax = p_alpha(2) + 20;
x = [ltf, htf, htf, ltf];
y = [0 0 ymax ymax];
g3 = patch(x, y, [114 114 114]./255,'FaceAlpha',.3, 'EdgeColor', 'none');
hold on
g1 = plot(fEO, pxxEO, 'LineWidth', 1.2);
g2 = plot(fEC, pxxEC, 'LineWidth', 1.2);
line([iaf iaf], [-5 ymax], 'Color', 'k'); % Plot Vertical Line - IAF
hold on
text(double([ltf iaf htf]), [ymax*0.5 ymax*0.7 ymax*0.9], {strcat('\leftarrow LTF = ',num2str(ltf)),...
    strcat('\leftarrow IAF = ',num2str(iaf)),strcat('\leftarrow HTF = ',num2str(htf))},'FontSize',8)

ylim([0 ymax])
xlim([1 18])
xlabel('Frequency (Hz)')
ylabel('Power (uV^2/Hz)')
hold off


