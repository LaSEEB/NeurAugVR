function EEG = epoch_and_merge(EEG, label)
events = EEG.event(strcmp({EEG.event(:).label}, label));
start = double([events(:).latency]);
finish = round(start + EEG.srate * [events(:).duration]);
EEG = pop_select(EEG, 'point', [start', finish']);
end