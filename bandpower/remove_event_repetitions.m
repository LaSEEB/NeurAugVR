function EEG = remove_event_repetitions(EEG)
% Obs.: An event with the same type and latency is considered repeated
types = {EEG.event(:).type};
latencies = [EEG.event(:).latency];
[~,ia,~] = unique(table(latencies(:), types(:)));
EEG.event = EEG.event(ia);
end