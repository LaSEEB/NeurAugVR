clc, clear, close all
%% Path
rfol = 'C:/Users/guta_/Desktop/Data Analysis/';             % Root folder
addpath(strcat(rfol,'Libraries/eeglab2021.0'))
varsbefore = who; eeglab; close; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})  % Run EEGLAB but remove variables and close figure

%% Parameters
fs = 250;  % Sampling frequency
fi = 10;  % Sinusoid frequency

%% Events
event = struct();
t = 0;
c = 1;
for i = 1:3  % Repetitions   
    [event,t,c] = add_event(event,fs,t,c,'Base',15);
    
    for j = 1:3  % Left trials
        [event,t,c] = add_event(event,fs,t,c,'Cross',5);
        [event,t,c] = add_event(event,fs,t,c,'L',5);
        [event,t,c] = add_event(event,fs,t,c,'EoT',1.25);
    end
    
    for k = 1:3  % Right trials
        [event,t,c] = add_event(event,fs,t,c,'Cross',5);
        [event,t,c] = add_event(event,fs,t,c,'R',5);
        [event,t,c] = add_event(event,fs,t,c,'EoT',1.25);
    end
end
[event,t,c] = add_event(event,fs,t,c,'EoE',6);

%% Data
times = (0:(t*fs-1))/fs;
data = sin(2*pi*fi*times);

% t = 0;
% [data,t] = modulate_sinusoid(data,fs,t,30,[100,100]);

t = 0;
for i = 1:3  % Repetitions 
    [data,t] = modulate_sinusoid(data,fs,t,15,[100,100]);
    
    for j = 1:6  % Left and right trials
        [data,t] = modulate_sinusoid(data,fs,t,5,[100,100]);
        [data,t] = modulate_sinusoid(data,fs,t,0.5,[100,20]);
        [data,t] = modulate_sinusoid(data,fs,t,4.5,[20,20]);
        [data,t] = modulate_sinusoid(data,fs,t,0.5,[20,100]);
        [data,t] = modulate_sinusoid(data,fs,t,0.75,[100,100]);
    end
end
[data,t] = modulate_sinusoid(data,fs,t,6,[100,100]);
% 
EEG = eeg_emptyset;
EEG.srate = fs;
EEG.event = event;
EEG.times = times*1000;
EEG.data = data;
EEG.chanlocs = struct('labels','C4');

%% Plot
figure
plot(EEG.times,EEG.data)
hold on
markers = {'Base','Cross','L','R','EoT','EoE'};
cols = jet(numel(markers));
mhs = gobjects(1,numel(markers));
ylim([-120,120])
for m = 1:numel(markers)
    mlats = [EEG.event(strcmp({EEG.event(:).type},markers{m})).latency];
    for mlat = mlats
        mh = plot(EEG.times(mlat)*[1,1],ylim,'Color',cols(m,:));
        hold on
    end
    mhs(m) = mh;
end
legend(mhs, markers, 'location', 'bestoutside')

% xlim([19700,20900])
% ylim([-120,120])
%% Save
% save('EEGsimu01.mat', 'EEG');
