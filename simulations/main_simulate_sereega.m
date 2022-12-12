clear, clc, close all
%% Add paths
rfol = 'C:/Users/guta_/Desktop/Data Analysis/';             % Root folder
cfol = fileparts(matlab.desktop.editor.getActiveFilename);  % Current file folder
restoredefaultpath; clear RESTOREDEFAULTPATH_EXECUTED;      % Remove any added paths
cd(cfol);                                                   % Change directory to the executing file's directory
addpath(genpath(strcat(rfol,'Libraries/SEREEGA/')))         % Add libraries. DOWNLOAD SEREEGA FROM: https://github.com/lrkrol/SEREEGA
addpath(strcat(rfol,'Libraries/eeglab2021.0'))
clear rfol cfol
varsbefore = who; eeglab; close; varsnew = []; varsnew = setdiff(who, varsbefore); clear(varsnew{:})  % Run EEGLAB but remove variables and close figure

%% Set
protocol = struct('n', 20,'srate',250,'length',8000,'prestim',3000,'marker',"S  7");  % Decided not to make 2x1 struct, but 2x1 'marker', since the EEGs will get merged!
leadfield = struct('montage', 'BioSemi32');
component = {struct('type','ersp','amplitude',3,'modMinRelAmplitude',0.5,'frequency',[8 10 10 12],'modulation','invburst','modLatency',5500,'modTaper',0.8,'modWidth',5000,'sources',{{[ 30, -17, 49]}});
             struct('type','noise','amplitude',3,'color','pink','amplitudeDv',1,'nsources',30,'spacing',20)
            };
sensornoise = 2;

%% Generate
EEGs = {};
for j = 1:numel(protocol.marker)
    %% Set channels, sources, and lead field
    lf = lf_generate_fromnyhead('montage', leadfield.('montage'));
    
    %% Make source orientation perpendicular to scalp
    perp_oris = utl_get_orientation_pseudoperpendicular(1:size(lf.pos,1),lf);
    lf.orientation = perp_oris;
    
    %% Components
    components = [];
    for k = 1:numel(component)
        comp = utl_check_class(component{k});
        switch comp.type
            case 'ersp', comp = ersp_component(lf, comp, comp.sources{j});
            case 'noise',comp = noise_component(lf, comp, comp.nsources(j), comp.spacing(j));
        end
        components = [components, comp];
    end
    
    %% Plot (optional)
%     ersp_source = lf_get_source_nearest(lf, component{1}.sources{j});  % Could perhaps automatically know it's component{1}
%     plot_source_location(ersp_source, lf, 'mode', '3d');
%     hold on
%     
%     p1 = [lf.pos(ersp_source,1), lf.pos(ersp_source,2), lf.pos(ersp_source,3)];
%     o1 = [lf.orientation(ersp_source,1), lf.orientation(ersp_source,2), lf.orientation(ersp_source,3)]*15;
%     mArrow3(p1,p1 + o1,'stemWidth', 0.4,'color',[1,0,0]);
%     hold on
%     
%     noise_nsources = component{2}.nsources;  % Could perhaps automatically know it's component{2}
%     noise_comp = components(2);
%     for k = 1:noise_nsources
%         s1 = [lf.pos(noise_comp.source(k),1), lf.pos(noise_comp.source(k),2), lf.pos(noise_comp.source(k),3)];
%         o1 = [lf.orientation(noise_comp.source(k),1), lf.orientation(noise_comp.source(k),2), lf.orientation(noise_comp.source(k),3)]*15;
%         scatter3(s1(1),s1(2),s1(3),5,'k','filled')
%         hold on
%         mArrow3(s1,s1+o1,'stemWidth', 0.4,  'tipWidth',1,'color',[0.4,0.4,0.4]);
%         hold on
%     end
    
    %% Data
    scalpdata = generate_scalpdata(components, lf, protocol,'sensornoise',sensornoise);

    %% EEG
    protocol_marker = protocol;
    protocol_marker.marker = protocol.marker(j);
    EEG = utl_create_eeglabdataset(scalpdata, protocol_marker, lf);
    
    %% Store
    EEGs{j} = EEG;
end

%% Extract one single EEG
if numel(EEGs) == 1
    EEG = EEGs{1};
elseif numel(EEGs) == 2  % This occurs when Left and Right trials were generated!
    EEG = merge_EEGs(EEGs);
else
    fprintf('More than 2 EEGs found!\n')
end
