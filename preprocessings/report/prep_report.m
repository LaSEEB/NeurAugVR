function prep_report(EEG)

set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'xtick',[])
set(gca,'xticklabel',[])

if ~isempty(EEG.preproc)
    fields = fieldnames(EEG.preproc);
    
    xti = [];
    xlab = {};
    for j = 1:numel(fields)
        ymax = size(EEG.data,1);
        if ismember('rej_chans',fields)
            ymax = ymax + numel(EEG.preproc.rej_chans);
        end
        switch fields{j}
            case 'interp_chans', [xti(end+1), xlab(end+1)] = plot_report_interp_chans(EEG, numel(xti)+1);  % interp_chans
            case 'rej_chans', [xti(end+1), xlab(end+1)] = plot_report_rej_chans(EEG, numel(xti)+1);  % interp_chans
            case 'rej_comps', [xti(end+1), xlab(end+1)] = plot_report_rej_comps(EEG, numel(xti)+1,ymax);
            case 'rej_trials', [xti(end+1), xlab(end+1)] = plot_report_rej_trials(EEG, numel(xti)+1,ymax);
            case 'rej_segments', [xti(end+1), xlab(end+1)] = plot_report_rej_segments(EEG, numel(xti)+1,ymax);
        end
    end
    xlim([0,xti(end)+1])
    ylim([0,ymax])
    
    xticks(xti)
    xticklabels(xlab)
else
    
end
end
