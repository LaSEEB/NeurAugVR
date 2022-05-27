function [xti, xlab] = plot_report_rej_comps(EEG, xi, ymax)
colors = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880];

compst = 1:size(EEG.icaweights,2);
compsr = EEG.preproc.rej_comps;
compsk = compst(~ismember(compst,compsr));

colids = [1*ones(1,numel(compsr)), 2*ones(1,numel(compsk))];
compssort = [sort(compsr), sort(compsk)];

n = numel(compssort);

xtemp = [0,1,1,0] -1 + xi*1.5;
% ytemp = [0,0,1,1];

if ymax
    ys = linspace(0,ymax, n + 1);
end

for i = 0:n-1
    x = xtemp;
%     y = ytemp + i;
    y = [ys(i+1), ys(i+1), ys(i+1+1), ys(i+1+1)];
    
    patch('XData',x,'YData',y,'FaceColor',colors(colids(i+1),:))
    hold on
    mx = mean(x(1:2));
    my = mean(y(2:3));
    % text(mx,my,chans{i+1},'HorizontalAlignment','center','FontSize',6)
    text(mx,my,sprintf('%d',compssort(i+1)),'HorizontalAlignment','center','FontSize',8,'FontUnits','normalized')
    hold on
end
% xlim([0,3])
% ylim([0,n])
na = numel(compssort);
ni = numel(compsr);
rep = sprintf('%d / %d (%d%%)',ni,na,round(ni/na*100));
% Option 1:
% set(gca,'xtick',3,'xticklabel',{sprintf('Rej. comps%s%s','\newline',rep)})

xti = mx;
xlab = {sprintf('Rej. comps%s%s','\newline',rep)};

end