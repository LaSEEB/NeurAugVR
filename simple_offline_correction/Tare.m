function EEG = Tare(EEG)

EEG.times = EEG.times - EEG.times(1);
EEG.xmin = EEG.times(1)/1000;
EEG.xmax = EEG.times(end)/1000;

end