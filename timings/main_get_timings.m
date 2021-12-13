%% Add
clc, clear, close all
addpath('C:/Users/guta_/Documents/MATLAB/eeglab2021.0')
varsbefore = who; eeglab; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})

%% Set
dfol = 'C:\Users\guta_\Desktop\Data Analysis\Functions\fMRI Analysis\folder1\';
sfol = 'C:\Users\guta_\Desktop\Data Analysis\Functions\fMRI Analysis\folder2\';
markers = {'S  1','S  2','S  5','S  7','S  8','S 10','S 11','S 12'};
scan_marker = 'Scan Start';
lim_markers = {'S  1', 'S 12'};

% For trial onsets
target_markers = {'S  7', 'S  8'};
condi_markers = {''};
condi = @(events, e, condi_marker) true;
end_markers = 'S 10';
snames = {'_trialTimingsL', '_trialTimingsR'};

% For block onsets
% target_markers = {'S  5'};
% condi_markers = {'S  7', 'S  8'};
% condi = @(events, e, condi_marker) strcmp(events(min(e+1, numel(events))).type, condi_marker);
% end_markers = {'S  2', 'S 11'};
% snames = {'_blockTimingsL', '_blockTimingsR'};

%% Find
cols = jet(numel(markers));  % Colors
file_struct = dir(strcat(dfol, '*.vhdr'));

for f = 1:numel(file_struct)
% for f = 1:1
    [~,dname,dext] = fileparts(file_struct(f).name);
    dnamex = strcat(dname, dext);
    EEG = pop_loadbv(dfol, dnamex, [], []);
    
    TR = 1.26;
    thr = 0.02;
    [~, EEGMARK, ~] = volume_markers(EEG, TR, thr);
    start = EEGMARK.times(EEGMARK.event(find(strcmp({EEGMARK.event(:).type}, scan_marker), 1, 'first')).latency);
    
    % Plot for debugging
    figure
    chn = 5;
    plot(EEG.times/1000, EEG.data(chn,:))
    hold on
    plot(EEG.times([EEGMARK.event(strcmp({EEGMARK.event(:).type}, scan_marker)).latency])/1000, EEG.data(chn, [EEGMARK.event(strcmp({EEGMARK.event(:).type}, scan_marker)).latency]), 'r*')
    hold on
    for m = 1:numel(markers)
        mlats = [EEG.event(strcmp({EEG.event(:).type},markers{m})).latency];
        for mlat = mlats
            plot(EEG.times(mlat)*[1,1]/1000,ylim, 'Color',cols(m,:), 'LineStyle', '-')
            hold on
        end
    end
    plot(EEG.times(EEG.event(find(strcmp({EEG.event(:).type}, lim_markers{1}), 1, 'first')).latency)/1000*[1,1], ylim, '--r')
    hold on
    
    events = EEG.event(ismember({EEG.event(:).type}, markers));
    % This is an extra-careful step to limit protocol event between the
    % first and last of one run (useful if the run has been repeated in the end, because of a Recview replay for example)
    lims = [find(strcmp({events(:).type}, lim_markers{1}), 1, 'first'), find(strcmp({events(:).type}, lim_markers{2}), 1, 'first')];
    events = events(lims(1):lims(2));
    
    for m = 1:max(numel(target_markers), numel(condi_markers))
        target_marker = target_markers{min(m, numel(target_markers))};
        condi_marker = condi_markers{min(m, numel(condi_markers))};
        
        mat = [];
        found = false;
        for e = 1:numel(events)
            if ~found && strcmp(events(e).type, target_marker) && condi(events, e, condi_marker)
                onset = EEG.times(events(e).latency) - start;
                found = true;
            end
            
            if found && any(strcmp(events(e).type, end_markers))
                dur = EEG.times(events(e).latency) - start - onset;
                mat = [mat; onset/1000, dur/1000, 0];
                found = false;
            end
            
        end
        % Plot for debugging
        if m == 1, style = '--k'; elseif m == 2, style = '-k'; end
        for i = 1:size(mat,1)
            plot(mat(i,1)*[1,1]+start/1000,ylim,style)
            hold on
        end
        writematrix(mat,sprintf('%s%s%s.txt',sfol, dname, snames{m}),'Delimiter','space')
    end
    title(sprintf('%s',dname), 'interpreter' ,'none')
end
