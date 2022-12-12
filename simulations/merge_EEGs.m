function EEG = merge_EEGs(EEGs)

EEG_Ltr = EEGs{1};
EEG_Rtr = EEGs{2};

%% Combine R/L EEGs into one
rand_vec = randperm(size(EEG_Ltr.data,3) + size(EEG_Rtr.data,3));
events = [EEG_Ltr.event, EEG_Rtr.event];
data = cat(3,EEG_Ltr.data, EEG_Rtr.data);
epochs = [EEG_Ltr.epoch, EEG_Rtr.epoch];

bias_lat = EEG_Ltr.event(end).latency;
bias_tim = EEG_Ltr.event(end).init_time;
bias_eve = EEG_Ltr.epoch(end).event;

events_temp = events;
data_temp = data;
epochs_temp = epochs;
for i = 1:size(rand_vec,2)
    events_temp(i).type = events(rand_vec(i)).type;
    data_temp(:,:,i) = data(:,:,rand_vec(i));
    epochs_temp(i).eventtype = epochs(rand_vec(i)).eventtype;
end
last_Ltr_lat = size(EEG_Ltr.data,2)*size(EEG_Ltr.data,3);
for i = size(EEG_Ltr.event,2)+1: size(events,2)
    events_temp(i).epoch = events_temp(i).epoch + bias_eve;
    events_temp(i).latency = events_temp(i).latency + last_Ltr_lat;
    events_temp(i).init_time = (events_temp(i).latency-1)/EEG_Ltr.srate;
    epochs_temp(i).event = events_temp(i).epoch;
    epochs_temp(i).eventinit_time = events_temp(i).init_time;
end

% Store combined data
EEG = EEG_Ltr;
EEG.data = data_temp;
EEG.event = events_temp;
EEG.epoch = epochs_temp;
EEG.trials = size(data_temp,3);

%% Deepoch
EEG = Deepoch(EEG);

end