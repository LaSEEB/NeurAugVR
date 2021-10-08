% ecgPeakDetection() - ECG peak detector based on the Tompkin's algorithm
%
% Usage:
%  >> [ ecg_at_peaks, Rpeaks ] = ecgPeakDetection(ecg_signal, fs)
%  
% Inputs:
%   ecg_signal    = 1D array 
%   fs            = data sampling rate (Hz)
%
% Outputs:
%    bpm    = beats per minute
%    Rpeaks = time-samples of the R peaks
%
% signal -> pode ser o ECG apenas, ou uma estrutura que tenha um campo "data" 
% e cuja ultima entrada seja o sinal de ECG (talvez seja mais claro olhando para o codigo na linha 17)
%
% fs -> frequencia de amostragem
%
% freq_band [ min_hz max_hz ] -> banda de frequencia para a qual queres que o 
% sinal seja filtrado (se vazio, [], entao nao faz filtragem ao ECG)
%
% reverse -> faz a detecao picos tambem do fim para o inicio, em vez de ser 
% %apenas do inicio para o fim (acrescentei esta parte porque reparava que o 
% %algoritmo falhava algumas vezes no inicio); depois, e feita uma combinacao das detecoes


function [ ecg, bpm, p, R_struct ] = ecgPeakDetection_v4(signal, fs, freq_band, reverse)

if isstruct(signal)
    ECG = signal.data(end, :);
else
    ECG = signal;
end

if ~isempty(freq_band)
    ECG = eegfilt(ECG, fs, freq_band(1), 0);
    ECG = eegfilt(ECG, fs, 0, freq_band(2));
end

if reverse
    ECG = vertcat(ECG, fliplr(ECG));
end   

for j = 1:size(ECG, 1)
    ecg = ECG(j, :);

    time = (1:length(ecg)) ./ fs;
    
    dp = 4; %4th order derivative
    dv = ecg(dp + 1:end) - ecg(1:end - dp);
    dv2 = dv .^ 2;
    
    sp = floor(0.100 * fs); %integration window of 100 ms
    
    sa = zeros(1, length(dv2));
    sb = zeros(1, length(dv2));
    osb = dv2(1) * 2; %initial threshold
    
    for i = 1:length(dv2)
        sa(i) = sum(dv2((i - min(sp, i)) + 1:i)) / sp; %falta resolver sa(1) = 0
        osb = osb + ((dv2(i) - osb) / 512.);
        sb(i) = osb;
    end
    
    sa = [0 sa(1:end - 1)];
    sb = sb .* 2; %ver no codigo do python o MULTIPLY
    
    c = sa > sb;
    c(1) = 0; c(end) = 0;
    
    dc = diff(c);
    ss = find(dc > 0); %the threshold is crossed
    se = find(dc < 0);
    
    % QRS Duration: between 40ms and 300ms
    
    mdt = 0.04 * fs;
    Mdt = 0.3 * fs;
    
    v = [];
    p = [];
    
    for i = 1:length(ss)
        si = ss(i); %start instant
        ei = se(i); %end instant
        
        [ ~, argmax ] = max(ecg(si:ei));
        ir = si + argmax;
        p = [p, ir];
        
        tr = time(ir);
        v = [v, tr];
    end
    
    %%%% change this part:
    p = qrscorrect(p, ecg, fs);
    
    thr04 = 0.4;
    c = find(diff(p) < (mean(diff(p)) * thr04));
    c_ = find(diff(p) < (mean(diff(p)) * thr04)) + 1;
    
    if isempty(c)
        p = p - 1;
        w = int64(mean(diff(p)) * 0.2);
        shift = zeros(1, length(p));
        
        for i = 1:length(p)
            if p(i) - w < 1
                w = length(1:p(i));
            elseif p(i) + w > length(ecg)
                w = length(p(i):length(ecg));
            end
            
            ecg_w = ecg((p(i) - w):(p(i) + w));
            [ ~, sh ] = max(ecg_w);
            shift(i) = sh - idivide(int64(length(ecg_w)), int64(2)) - 1;
        end
        
        p = p + shift;
        
        p(diff(p) == 0) = [];
    end
    
    while ~isempty(c)
        indexes = [];
        for i = 1:length(c)
            if ecg(p(c(i))) <= ecg(p(c_(i)))
                indexes = [indexes c(i)];
            else
                indexes = [indexes c_(i)];
            end
        end
        % v(indexes) = [];
        p(indexes) = [];
        
        p = p - 1;
        w = int64(mean(diff(p)) * 0.2);%0.2
        shift = zeros(1, length(p));
        
        for i = 1:length(p)
            if p(i) - w < 1
                w = length(1:p(i));
            elseif p(i) + w > length(ecg)
                w = length(p(i):length(ecg));
            end
            ecg_w = ecg((p(i) - w):(p(i) + w));
            [ ~, sh ] = max(ecg_w);
            shift(i) = sh - idivide(int64(length(ecg_w)), int64(2)) - 1;
        end
        
        p = p + shift;
        p(diff(p) == 0) = [];
        
        c = find(diff(p) < (mean(diff(p)) * 0.6));%0.6
        c_ = find(diff(p) < (mean(diff(p)) * 0.6)) + 1;%0.6
    end
    
    P{j} = p;
end


%%%%% change until here

if reverse
    p1 = P{1}; p2 = P{2};
    l = length(ECG):-1:1; p2 = fliplr(l(p2));
    p = [ p2(p2 < p1(1)), p1 ];
else
    p = P{1};
end

for a = 1:length(p)
    R_struct(a).latency  = p(a);
    R_struct(a).duration = 1;
    R_struct(a).channel  = 0;
    R_struct(a).type     = 'RT';
    R_struct(a).bvtime   = [];
    R_struct(a).bvmknum  = a;
    R_struct(a).urevent  = a;
    R_struct(a).code     = 'Comment';
end

time_sec = length(ECG) / fs;
time_min = time_sec / 60;
bpm = length(p) / time_min;
ecg = ECG(1, :);
