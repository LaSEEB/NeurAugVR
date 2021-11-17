function plot_markers(EEG, markers, chn)

cols = jet(numel(markers));
for m = 1:numel(markers)
    lats = [EEG.event(strcmp({EEG.event(:).type},markers{m})).latency];
    
    if isempty(chn)  % If no chan was specified, draw vertical lines
        for lat = lats
            plot(EEG.times(lat)*[1,1]/1000,ylim,'Color',cols(m,:), 'LineStyle', '-')
            hold on
        end
        
    elseif isnumeric(chn) % If a chan was specified, draw marks
        plot(EEG.times(lats)/1000, EEG.data(chn,lats), 'Color', 'r', 'LineStyle', 'none', 'Marker', '*');
        hold on
    end
    
end
end