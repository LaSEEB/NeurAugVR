function EEG = add_event(EEG,type,lat)

% Make events vertical
EEG.event = EEG.event(:);

N = numel(EEG.event);

EEG.event(N+1).type = type;
EEG.event(N+1).latency = lat;

[~,I] = sort([EEG.event(:).latency]);
EEG.event = EEG.event(I);
for j = 1:size(EEG.event,1)
    EEG.event(j).urevent = j;
end

end