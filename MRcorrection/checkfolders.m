function checkfolders(path)

markers = [path, '/Markers'];
dws = [path, '/DWS'];
qrs = [path, '/QRS'];
ica = [path, '/ICA'];
PAfiltered = [path, '/PAfiltered'];
powerfreq = [path, '/PowerFreq'];
ratioC = [path, '/ratioC'];
OriginalMat = [path, 'Original_mat'];
Figures = [path, 'Figures'];

%Figures

 if ~exist(Figures, 'dir')
        mkdir(Figures)
 end


%Orginal Mat

    if ~exist(OriginalMat, 'dir')
        mkdir(OriginalMat)
    end

%markers
    if ~exist(markers, 'dir')
       mkdir(markers)
    end
    
%GA filtered, cut and downsampled    
    if ~exist(dws, 'dir')
       mkdir(dws)
    end
    
% QRS - R markers    
    if ~exist(qrs, 'dir')
       mkdir(qrs)
    end
    
% ICA needed for PROJIC    
    if ~exist(ica, 'dir')
       mkdir(ica)
    end
    
% powerfreq
    if ~exist(powerfreq, 'dir')
       mkdir(powerfreq)
    end
    
% ratioC
    if ~exist(ratioC, 'dir')
       mkdir(ratioC)
    end
    
% PAfiltered
    if ~exist(PAfiltered, 'dir')
       mkdir(PAfiltered)
    end