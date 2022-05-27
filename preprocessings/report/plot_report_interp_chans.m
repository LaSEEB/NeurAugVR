function [xti, xlab] = plot_report_interp_chans(EEG, xi)
colors = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880];
chans = {EEG.chanlocs.labels};
% chans = [{EEG.chanlocs.labels}, {EEG.preproc.rejchans}];
chansi = EEG.preproc.interp_chans; % {EEG.preproc.intchans};

if isfield(EEG.preproc,'rej_chans')
    chans = [chans, EEG.preproc.rej_chans];
end

chansni = [chans(~ismember(chans,chansi)),];
chanssort = [flip(sort(chansi)), flip(sort(chansni))];
chancolids = [1*ones(1,numel(chansi)), 2*ones(1,numel(chansni))];
n = numel(chanssort);

xtemp = [0,1,1,0] -1 + xi*1.5;
ytemp = [0,0,1,1];
for i = 0:n-1
    x = xtemp;
    y = ytemp + i;
    patch('XData',x,'YData',y,'FaceColor',colors(chancolids(i+1),:))
    hold on
    mx = mean(x(1:2));
    my = mean(y(2:3));
    % text(mx,my,chans{i+1},'HorizontalAlignment','center','FontSize',6)
    text(mx,my,chanssort{i+1},'HorizontalAlignment','center','FontSize',8,'FontUnits','normalized')
    hold on
end
% xlim([0,3])
% ylim([0,n])
na = numel(chanssort);
ni = numel(chansi);
rep = sprintf('%d / %d (%d%%)',ni,na,round(ni/na*100));
% Option 1:
% set(gca,'xtick',mx,'xticklabel',{sprintf('Interp. chans%s%s','\newline',rep)})

xti = mx;
xlab = {sprintf('Interp. chans%s%s','\newline',rep)};
end