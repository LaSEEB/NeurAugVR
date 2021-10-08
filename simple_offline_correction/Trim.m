function EEGCUT = Trim(EEGX, lims, units)

if strcmp(units, 'time')
    [idx1, idx2] = Get_ROI_indexes(EEGX, lims);
elseif strcmp(units, 'idx')
    idx1 = lims(1);
    idx2 = lims(2);
elseif strcmp(units, 'marker')
    idx1 = EEGX.event(find(strcmp({EEGX.event(:).type},lims(1)),1,'first')).latency;
    idx2 = EEGX.event(find(strcmp({EEGX.event(:).type},lims(2)),1,'last')).latency - 1000; % -1000, CAREFUL
end

% cut data
EEGCUT = EEGX;
EEGCUT.data = EEGX.data(:, idx1:idx2);
EEGCUT.times = EEGX.times(:, idx1:idx2);
EEGCUT.xmin = EEGCUT.times(1)/1000;
EEGCUT.xmax = EEGCUT.times(end)/1000;
% EEGCUT.time = 1:length(EEGCUT.data)/EEGGAR.srate;
% EEGCUT.times = (0:length(EEGCUT.data)-1)/EEGX.srate;
EEGCUT.pnts  = size(EEGCUT.data,2);

if ~isempty(EEGX.event)
    event_idx_mask = [EEGX.event(:).latency] >= idx1 & [EEGX.event(:).latency] <= idx2;
    % correct events structure taking this change into account 
    EEGCUT.event = EEGX.event(event_idx_mask);
%     new_latencies = [EEGCUT.event(:).latency] - EEGCUT.event(1).latency + 1;
    new_latencies = [EEGCUT.event(:).latency] - idx1 + 1;
    new_latencies_cell = num2cell(new_latencies);
    [EEGCUT.event.latency] = deal(new_latencies_cell{:});
    
    for i = 1:size(new_latencies, 2)
        EEGCUT.event(i).urevent = i;
    end
    
end

%% To visualize:
% chan = 32;
% figure;
% plot(EEGX.times/1000, EEGX.data(chan,:))
% hold on
% plot(EEGCUT.times/1000, EEGCUT.data(chan,:))
% % titleplot = strcat({['Artifact Deletion', ' - ', subject],['Channel: ', EEGX.chanlocs(chan).labels]});
% % title(titleplot)
% legend('Original EEG', 'Cut EEG')
% ylabel(' Amplitude [uV] ')
% xlabel('Time [s]')
% hold off

end