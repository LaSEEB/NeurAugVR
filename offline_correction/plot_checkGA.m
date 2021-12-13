function plot_checkGA(filename, EEG, EEGGAR, TR, SMS, nr_slices, pnts, chans, locfilename)
% plot_checkGA(filename, EEG, EEGGAR, TR, SMS, nr_slices, pnts, chans, locfilename)
%
% Creates figures to evaluate GA correction quality (it does not save them!)
%
% filename = name of the file, for plotting
% plot_checkGA(filename, EEG, EEGGAR, TR, SMS, nr_slices, pnts, chans, locfilename)
% EEG = uncorrected EEG with volume markers specified as EEG.volmarkers
% EEGGAR = EEG corrected for the GA artifact
% TR = repetition time in seconds - 1.26 s
% SMS = simultaneous multi-slice acquisition factor - 3;
% nr_slices = total number of slices - 60;
% pnts = index of the first and last point of the acquisition - [first_vol last_vol+TR*EEGGAR.srate-1];
% chans = channels of interest for the plots - [5, 6, 18]; %C3, C4, Cz
% locfilename = file with electrode 

%% Distribution of amplitude values before and after GA correction for all channels
figure('pos', [50 50 1100 600]);
h = suptitle({filename, '  ', '  '});
h.Interpreter = 'none';
subplot(2,1,1)
boxplot([EEG.data'])
set(gca, 'XTickLabels', {EEG.chanlocs.labels}, 'XTickLabelRotation', 45)
ylabel('Amplitude [uV]')
title('Before GA correction')

subplot(2,1,2)
boxplot([EEGGAR.data(:, pnts(1):pnts(end))'])
set(gca, 'XTickLabels', {EEG.chanlocs.labels}, 'XTickLabelRotation', 45)
ylabel('Amplitude [uV]')
title('After GA correction')

%% Time and frequency analysis for the channels specified in "chans"

t = (0:1:EEG.pnts-1)./EEG.srate;
t_vol1 = pnts(1)/EEG.srate;
t_volend = pnts(2)/EEG.srate; 

f_GA = 1/(TR/(nr_slices/SMS));
L = length(EEG.data(:, pnts(1):pnts(end)));
freq = (0:L-1)*(EEG.srate/L);
[~, f_40] = min(abs(freq-40)); % find the index that corresponds to 40 Hz
[~, f_01] = min(abs(freq-0.10)); % find the index that corresponds to 0.1 Hz
[~, f_GA_l] = min(abs(freq-(f_GA-5))); % find the index that corresponds to 40 Hz
[~, f_GA_r] = min(abs(freq-(f_GA+5))); % find the index that corresponds to 0.1 Hz


for c = 1:length(chans)
    chan = chans(c);
    figure('pos', [50 50 1100 600]);
    h = suptitle([filename, '; channel: ', EEG.chanlocs(chan).labels, '; f_GA = ', num2str(round(f_GA,2)), ' Hz']);
    h.Interpreter = 'none';
    
    % Whole signal: before vs after GA correction
    subplot(3,2,[1 2])
    plot(t, EEG.data(chan, :))
    hold on
    plot(t, EEGGAR.data(chan, :))
    line([pnts(1)/EEG.srate pnts(1)/EEG.srate], ylim, 'LineStyle', '--', 'Color', 'k')
    line([pnts(2)/EEG.srate pnts(2)/EEG.srate], ylim, 'LineStyle', '--', 'Color', 'k')
    ylabel('Amplitude [uV]')
    xlabel('Time [s]')
    legend('Original', 'GAremoved')
    
    % Power spectrum 
    pEEG = abs(fft(EEG.data(chan, pnts(1):pnts(end)))).^2/L;
    pEEGGAR = abs(fft(EEGGAR.data(chan, pnts(1):pnts(end)))).^2/L;
    
    %%%before vs after GA correction - 0.1 to 40 Hz
    ymin = min(pEEGGAR(f_01:f_40));
    ymax = max(pEEGGAR(f_01:f_40));
    subplot(3,2,3)
    plot(freq, pEEG)
    hold on 
    plot(freq,pEEGGAR)
    line([f_GA f_GA], [ymin ymax], 'LineStyle', '--', 'Color', 'k')
    xlim([0.1 40])
    ylim([ymin ymax])
    xlabel('Frequency [Hz]')
    ylabel('Power [uV²]')
    legend('Original', 'GAremoved')
    
    %%% after GA correction - around f_GA
    ymin = min(pEEGGAR(f_GA_l:f_GA_r));
    ymax = max(pEEGGAR(f_GA_l:f_GA_r));
    subplot(3,2,4)
    plot(freq,pEEGGAR)
    
    %hold on;
    line([f_GA f_GA], [ymin ymax], 'LineStyle', '--', 'Color', 'k')
    xlim([f_GA-5 f_GA+5])
    ylim([ymin ymax])
    xlabel('Frequency [Hz]')
    ylabel('Power [uV²]')
    
    
    % Transition: before vs after GA correction
    tbefore = 15;
    tafter = 10;
    subplot(3,2,5)
    plot(t, EEG.data(chan, :))
    hold on
    plot(t, EEGGAR.data(chan, :))
    ylabel('Amplitude [uV]')
    xlabel('Time [s]')
    legend('Original', 'GAremoved')
    xlim([t_vol1-tbefore t_vol1+tafter])
    
    % After GA correction
    subplot(3,2,6)
    plot(t, EEGGAR.data(chan, :))
    ylabel('Amplitude [uV]')
    xlabel('Time [s]')
    xlim([t_vol1 t_vol1+tafter])
    
end

%% GA Artifact 


%%
nr_channels = 32;
slice_plot = 2;
slice_duration = TR/(nr_slices/SMS);
dif = zeros(nr_channels, slice_plot*slice_duration*EEG.srate);
vol = EEG.volmarkers; 

EEGepoch = pop_epoch(EEG, {  'Scan Start'  }, [0        1.26], 'newname', ...
    [filename, '_epochs'], 'epochinfo', 'yes');

% All channels

figure('pos', [50 50 1300 600]);
h = suptitle({[filename, ' - GA artifact - 2 slices - All channels'], '  ', '  '});
h.Interpreter = 'none';
for k = 1:nr_channels

    subplot(4, 8, k)
    data = reshape(EEGepoch.data(k,1:slice_plot*slice_duration*EEG.srate,:), slice_plot*slice_duration*EEG.srate, vol);
    plot(0:1/EEG.srate:slice_plot*slice_duration-1/EEG.srate, data, 'k')
    m = mean(data, 2);
    
    hold on
    plot(0:1/EEG.srate:slice_plot*slice_duration-1/EEG.srate, m, 'Linewidth', 1, 'Color', 'm')
    xlim([0 slice_plot*slice_duration])
    if ismember(k,[25, 26, 27, 28, 29, 30, 31, 32])
        xlabel('Time [s]')
    end
    if ismember(k, [1, 9, 17, 25])
        ylabel('Amplitude [uV]')
    end
    title([EEGepoch.chanlocs(k).labels])
    
    d = (data-repmat(m(1:slice_plot*slice_duration*EEG.srate), 1, vol)).^2;
    dif(k, :) = sum(d,2);

end

%% Differences
figure('pos', [50 50 1300 600]);
st = ':';
for n = 1:31
    if n > 10 & n < 20
        st = '-.';
    elseif n> 20
        st = '--';
    end       
    plot(0:1/EEG.srate:slice_plot*slice_duration-1/EEG.srate, dif(n, :)', 'Linewidth', 1, 'Linestyle', st)
    hold on
end
xlim([0 slice_plot*slice_duration])
xlabel('Time [s]')
ylabel('Sum (Trial-M)^2')
legend(EEG.chanlocs.labels)
title(['Trial difference - ', filename], 'Interpreter', 'none')

dif_norm = dif./(m'.^2);
sumdif = sum(dif_norm,2);

figure('pos', [50 50 1200 500]);
bar(1:31, sumdif(1:31))
set(gca, 'XTick', 1:31, 'XTickLabels', {EEG.chanlocs(1:31).labels}')
xlim([0 32])
ylabel('Sum along volumes and time points ((Trial-M)^2)/M^2')
title(['Trial difference sum along time - ', filename], 'Interpreter', 'none')

figure('pos', [50 50 1000 600]);
topoplot(sumdif, locfilename, 'electrodes', 'labels');
title(['Trial difference sum along time - ', filename], 'Interpreter', 'none')