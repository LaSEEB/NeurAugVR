function [xti, xlab] = plot_report_rej_segments(EEG, xi, ymax)

colors = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880];

segp = EEG.preproc.rej_segments;  % percent


xtemp = [0,1,1,0] -1 + xi*1.5;
x = xtemp;
y = [0, 0, ymax * segp / 100, ymax * segp / 100];
patch('XData',x,'YData',y,'FaceColor',colors(1,:))
hold on
y = [ymax * segp / 100, ymax * segp / 100, ymax, ymax];
patch('XData',x,'YData',y,'FaceColor',colors(2,:))
hold on

mx = mean(x(1:2));
xti = mx;
rep = sprintf('%d%%',round(segp));
xlab = {sprintf('Rej. segments%s%s','\newline',rep)};

end