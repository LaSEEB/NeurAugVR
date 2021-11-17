function EEG = undo_epochs(EEG)

% Reshape data
EEG.data = reshape(EEG.data, size(EEG.data,1), size(EEG.data,2)*size(EEG.data,3));

% To put 'boundary' events representing cuts, but commented since it gives
% an error if epoched again
% temp_array = EEG.event(1);
% for i = 1:EEG.trials-1
%     temp_array(i).type = 'boundary';
%     temp_array(i).latency = i*EEG.pnts-0.5;
%     temp_array(i).duration = NaN;
% end
% 
% EEG.event = [EEG.event, temp_array];
% [~,I] = sort([EEG.event(:).latency]);
% EEG.event = EEG.event(I);
% for i = 1:size(EEG.event,2)
%     EEG.event(i).urevent = i;
% end

% Fill rest of details
EEG.event = rmfield(EEG.event,{'epoch'});
EEG.urevent = rmfield(EEG.event,{'urevent'});
EEG.times = (0:(EEG.pnts*EEG.trials-1))/EEG.srate*1000;
EEG.xmin = EEG.times(1)/1000;
EEG.xmax = EEG.times(end)/1000;
EEG.pnts = size(EEG.data,2);
EEG.epoch = [];
EEG.trials = 1;

EEG = eeg_checkset(EEG);
end
