function plot_time(EEG, chan, extra_marker)
blue_col = [0 0.4470 0.7410];
signal = EEG.data(chan,:);
plot(EEG.times/1000, signal, 'Color', blue_col);
hold on
markers = {'S  1','S  2','S  5','L','R','S 10','S 11','S 12'};
cols = jet(numel(markers));
for m = 1:numel(markers)
    lats = [EEG.event(strcmp({EEG.event(:).type},markers{m})).latency];
    for lat = lats
        plot(EEG.times(lat)*[1,1]/1000,ylim,'Color',cols(m,:))
        hold on
    end
end
if ~isempty(extra_marker)
    vol_lats = [EEG.event(strcmp(extra_marker, {EEG.event(:).type})).latency];
    plot(EEG.times(vol_lats)/1000, signal(vol_lats), 'r*');
    hold on
end

end