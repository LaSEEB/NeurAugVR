function [fEO, fEC, pxxEO, pxxEC, f_alpha, p_alpha, band_alpha] = iab_features(data_EO, data_EC, fs, method)
%
% [fEO, fEC, pxxEO, pxxEC, f_alpha, p_alpha, band_alpha] = iab_features(data_EO, data_EC, fs, method)
%
% INPUT
%
% > data_EO - 1xm EEG signal, with m timepoints, recorded for eyes open
% > data_EC - 1xm EEG signal, with m timepoints, recorded for eyes closed
% > srate - sampling frequency
% > method - method to define the individual alpha band. Options:
% -- standard (LTF = 8 Hz, IAF = 10 Hz, HTF = 12 Hz);
% -- fbfw - fixed band fixed width (LTF = 6 Hz, IAF = 8 Hz, HTF = 12 Hz);
% -- ibfw - individual band fixed width
% -- ibiw - individual band individual width
% -- jsousa - Joana Sousa's method to adjust LTF and HTF according to IAF 
% position 
% -- crossings - nearest intersection points surrounding the IAF
% -- manual - interactive manual selection 
% 
% OUTPUT
%
% > fEO - frequency array for eyes open
% > fEC - frequency array for eyes closed
% > pxxEO - power spectral density (PSD) for eyes open using Welch's method
% > pxxEC - power spectral density (PSD) for eyes closed using Welch's method
% (window = 1024, overlap = 50%)
% > f_alpha - [ltf; iaf; htf] - frequencies in Hz;
% > p_alpha - [p_ltf; p_iaf; p_htf] - PSD at f_alpha frequencies 
% > band_alpha - [sump_EO; sump_EC] - PSD sum for eyes open and eyes closed



% Define pwelch's inputs
window = 1024; % segment length 
noverlap = window*0.5; % segment overlap 
nfft = 2^nextpow2(length(data_EO)); % number of DFT points

[pxxEO, fEO] = pwelch(data_EO, window, noverlap, nfft, fs);
[pxxEC, fEC] = pwelch(data_EC, window, noverlap, nfft, fs);

% find frequencies's indexes (alpha and whole eeg boundaries)
[~,f30] = min(abs(fEO-30));
[~,f8] = min(abs(fEO-8));
[~,f13] = min(abs(fEO-13));

% Select only the frequency range from 0.4 to 30 Hz
EO = [(fEO(1:f30))';((pxxEO(1:f30))')];
EC = [(fEC(1:f30))';((pxxEC(1:f30))')];


% Calculate the intersections between the EO and the EC curves
P = InterX(EO,EC);

% Determine the Individual Alpha Frequency (IAF), the peak within the
% 8-13Hz range in the EC curve
[iaf, p2] = findPeak(fEC,pxxEC,[8,12]);

iaf = double(iaf);


switch method
    case 'standard'
        iaf = 10;
        ltf = 8;
        htf = 12;
        
    case 'fbfw'
        iaf = 8;
        ltf = 6;
        htf = 12;
    
    case 'ibfw'
        ltf = iaf - 4;
        htf = iaf + 2;
        
    case 'ibiw'
        ltf = iaf * 0.6;
        htf = iaf * 1.2;
        
    case  'jsousa' 
        if iaf <= 10

            ltf = iaf-(1-abs(iaf-10)/10)*2;
            htf = ltf + 4;

        else

            htf = iaf+(1-abs(iaf-10)/10)*2;
            ltf = htf - 4;
        end
        
    case 'crossings'
        % Determine the high and low transition frequencies (nearest intersection points surrounding the IAF)
        inter = P(1,:);
        dist = inter-iaf;

        Left = max (dist(dist<0));
        Rigth = min (dist(dist>0));

        LTF_ind = find(dist==Left);
        HTF_ind = find(dist==Rigth);

        ltf = double(inter(LTF_ind));
        htf = double(inter(HTF_ind));

    case 'manual'
        figure
        plot(EO(1,:),EO(2,:),EC(1,:),EC(2,:),P(1,:),P(2,:),'ro')
        xlabel('Frequency (Hz)')
        ylabel('Power (uV²/Hz)')
        ylim([-5 300])
        legend('EO','EC')
        hold on
        
        disp('Select LTF and press Enter')
        [ltf, ~] = ginput;
        text(double(ltf), 10, {strcat('\leftarrow LTF = ',num2str(ltf))})
        line([ltf ltf], [-5 350], 'Color', 'k') % Plot Vertical Line - LTF
        
        disp('Select IAF and press Enter')
        [iaf, ~] = ginput;
        text(double([ltf iaf]), [10 50], {strcat('\leftarrow LTF = ',num2str(ltf)),...
        strcat('\leftarrow IAF = ',num2str(iaf))})
        line([iaf iaf], [-5 350], 'Color', 'c') % Plot Vertical Line - IAF
        
            
%         disp('Select IAF and press Enter')
%         [iaf1, ~] = ginput;
%         text(double([ltf iaf1]), [10 50], {strcat('\leftarrow LTF = ',num2str(ltf)),...
%         strcat('\leftarrow IAF - peak1 = ',num2str(iaf1))})
%         line([iaf1 iaf1], [-5 350], 'Color', 'c') % Plot Vertical Line - IAF
%         
%         disp('Select IAF and press Enter')
%         [iaf2, ~] = ginput;
%         text(double([ltf iaf2]), [10 50], {strcat('\leftarrow LTF = ',num2str(ltf)),...
%         strcat('\leftarrow IAF - peak2 = ',num2str(iaf2))})
%         line([iaf2 iaf2], [-5 350], 'Color', 'c') % Plot Vertical Line - IAF
        
        disp('Select HTF and press Enter')
        [htf, ~ ] = ginput;
        text(double([ltf iaf htf]), [10 50 30], {strcat('\leftarrow LTF = ',num2str(ltf)),...
        strcat('\leftarrow IAF = ',num2str(iaf)),strcat('\leftarrow HTF = ',num2str(htf))})
        %text(double([htf]), [30], {strcat('\leftarrow HTF = ',num2str(htf))})
        line([htf htf], [-5 350], 'Color', 'g') % Plot Vertical Line - HTF

        pause 
        disp('Press Enter to close the figure')
        close
end

[~,f_iaf] = min(abs(fEC-iaf));
[~,f_ltf] = min(abs(fEC-ltf));
[~,f_htf] = min(abs(fEC-htf));
p_iaf = double(pxxEC(f_iaf));
p_ltf = double(pxxEC(f_ltf));
p_htf = double(pxxEC(f_htf));

sump_EO = double(sum(pxxEO(f_ltf:f_htf)));
sump_EC = double(sum(pxxEC(f_ltf:f_htf)));
        
f_alpha = [ltf; iaf; htf];
p_alpha = [p_ltf; p_iaf; p_htf];
band_alpha = [sump_EO; sump_EC];
        