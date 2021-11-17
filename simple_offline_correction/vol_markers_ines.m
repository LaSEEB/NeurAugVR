function EEG  = vol_markers_ines(EEG, TR, thr, etype)

% highPassFilter = 1;  % Hz
% data = pop_eegfiltnew(EEG, highPassFilter, []);

%% 
data = EEG;
new_data = data;
fs = EEG.srate;
grad_thr = 200; % gradient trigger 200 uV/ms
nr_slices = 20; % 20 slices per volume
TRslice = (TR/nr_slices)*EEG.srate; % number of points of a slice
w = (1/EEG.srate)*1000; % time between two points in ms

% use average of all channels (except ECG, the last channel):
n_chan = data.nbchan; 
l_chan = 1:n_chan;
template = ismember(l_chan, 1:n_chan-1); 
av_chan = (template * [data.data]  /sum(template));

grad_ids = find(abs(gradient(av_chan, w)) > grad_thr); % selects only the points
                                                  % above gradient threshold

% figure;
% plot(EEG.times/1000, av_chan)
% hold on
% plot(EEG.times(grad_ids)/1000, ones(1, length(grad_ids)), '*')
% xlabel('Time [s]')
% ylabel('Amplitude [uV]')


diff_ids = diff(grad_ids); % distance in points between the selected indexes
val = unique(diff_ids);
cnt = histc(diff_ids, val);
counting = [val; cnt];


thr_diff = diff_ids(diff_ids <thr*fs); % distance between consecutive
                                       % indexes has to be lower than
                                       % TRslice. In this case, a
                                       % stricter value, TRslice/3 was
                                       % used

                                           
m = max(thr_diff);
ind = unique(thr_diff(thr_diff >= m - 10)); % from those, select only the 
                                           % ones for which the distance 
                                           % does not differ from the
                                           % maximum more than 10 points
                                        
x = ismember(diff_ids, ind); 

ids = grad_ids(find(x)+1); 

% GUS: tirei à força o último marker, parecia estar a dar erros
% % take out last point if there is no artifact after it, by comparing the
% % maximum signal amplitudes in a window after the marker with the period 
% % before the marker
% try 
%     
% last = ids(end);
% w_before = last-TR/2*EEG.srate:last-1; % window before last point
% w_after = last+(TR/2)*EEG.srate+1:last+TR*EEG.srate;% window TR/2 after last point
%                                                     %(the window does not start 
%                                                     %right after the last marker
%                                                     %because there might be brief 
%                                                     %large amplitude oscillations)
% 
% if max(data.data(1, w_after)) < 0.5 * max(data.data(1, w_before))
%     ids = ids(1:end-1);
% end
% 
% catch
% end

% ids = ids(1:end-1);

% CONFIRM TIMINGS!
timings = diff(ids); 

% exclude points that do not differ approximately TR from their neighbour
if ~all(abs(timings-TR*fs) < 100) 
    ids = [ids(abs(timings-TR*fs) < 100), ids(end)];
end

% Exclude dummy ids
ids = ids(3:end);

% Update EEG
temp_array = EEG.event(1);
for i = 1:size(ids,2)
    temp_array(i).type = etype;
    temp_array(i).latency = ids(i);
end
EEG.event = [EEG.event, temp_array];
[~,I] = sort([EEG.event(:).latency]);
EEG.event = EEG.event(I);
for i = 1:size(EEG.event,2)
    EEG.event(i).urevent = i;
end

% % Trim
% % trim_lims = [ids(3)+1,ids(end-1)-1]; % Trim from 3rd vol-marker to 2nd-to-last
% % trim_lims = [ids(3),ids(end-1)]; % Trim from 3rd vol-marker to 2nd-to-last (inclusive!)
% % trim_lims = [ids(3),size(EEG.data,2)]; % Trim from 3rd vol-marker to 2nd-to-last (inclusive!)
% trim_lims = [ids(3),ids(end)]; % Trim from 3rd vol-marker to 2nd-to-last (inclusive!)
% 
% % scan_ids = [ids(3),ids(end)];
% 
% 
% EEGtrim = Trim(EEG, trim_lims, 'idx');
% EEGtrim_begining = Trim(EEG, [1, trim_lims(1)-1], 'idx');
% EEGtrim_ending = Trim(EEG, [trim_lims(2)+1, size(EEG.data,2)], 'idx');
% 
% %% Final plot with volumes
% fig = 0;
% % fig = figure;
% % Deleted begining
% ids_begining = [EEGtrim_begining.event(strcmp(etype, {EEGtrim_begining.event(:).type})).latency];
% av_chan_begining = (template * [EEGtrim_begining.data]  /sum(template));
% plot(EEGtrim_begining.times/1000, av_chan_begining,'Color', [0.5,0.5,0.5])
% hold on
% plot(EEGtrim_begining.times(ids_begining)/1000, av_chan_begining(ids_begining), 'ko')
% hold on
% 
% % Preserved part
% blue_col = [0 0.4470 0.7410];
% ids = [EEGtrim.event(strcmp(etype, {EEGtrim.event(:).type})).latency];
% av_chan = (template * [EEGtrim.data]  /sum(template));
% ph1 = plot(EEGtrim.times/1000, av_chan, 'Color', blue_col);
% hold on
% 
% %% Plot_cues(EEGtrim)
% markers = {'S  1','S  2','S  5','L','R','S 10','S 11','S 12'};
% cols = jet(numel(markers));
% for m = 1:numel(markers)
%     lats = [EEGtrim.event(strcmp({EEGtrim.event(:).type},markers{m})).latency];
%     for lat = lats
%         plot(EEGtrim.times(lat)*[1,1]/1000,ylim,'Color',cols(m,:))
%         hold on
%     end
% end
% 
% % Deleted ending
% ids_ending = [EEGtrim_ending.event(strcmp(etype, {EEGtrim_ending.event(:).type})).latency];
% av_chan_ending = (template * [EEGtrim_ending.data]  /sum(template));
% plot(EEGtrim_ending.times/1000, av_chan_ending,'Color', [0.5,0.5,0.5])
% hold on
% ph3 = plot(EEGtrim_ending.times(ids_ending)/1000, av_chan_ending(ids_ending), 'ko');
% 
% % Volumes
% ph2 = plot(EEGtrim.times(ids)/1000, av_chan(ids), 'r*');
% hold on
% 
% xlabel('Time [s]')
% ylabel('Amplitude [uV]')
% title(['Number of volumes: ', num2str(length(ids)+1)])
% % legend([ph1,ph2,ph3],'Average EEG', 'Volume markers', 'Excluded volume markers')
% 
% %% Tare
% EEGtrim = Tare(EEGtrim);

end