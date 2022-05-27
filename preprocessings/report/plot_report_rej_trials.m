function [xti, xlab] = plot_report_rej_trials(EEG, xi, ymax)

trialnames = fieldnames(EEG.preproc.rej_trials);
nnames = numel(trialnames);

colors = lines(nnames+1);
colors(2,:) = [];
colors = [0.8500 0.3250 0.0980; colors];

ntvec = zeros(1,nnames);
nrvec = zeros(1,nnames);

nt = 0;

totaltrials = 0;
for i = 1:nnames
    trialsk = EEG.preproc.rej_trials.(trialnames{i}).kept;
    trialsr = EEG.preproc.rej_trials.(trialnames{i}).rej;
    colids = [1*ones(1,numel(trialsr)), (1+i)*ones(1,numel(trialsk))];
    trialssort = [sort(trialsr), sort(trialsk)];
    totaltrials = totaltrials + numel(trialssort);
end
if ymax
    ys = linspace(0,ymax, totaltrials + 1);
end

ttr = 0;
for i = 1:nnames
    
    trialsk = EEG.preproc.rej_trials.(trialnames{i}).kept;
    trialsr = EEG.preproc.rej_trials.(trialnames{i}).rej;
    colids = [1*ones(1,numel(trialsr)), (1+i)*ones(1,numel(trialsk))];
    trialssort = [sort(trialsr), sort(trialsk)];
    ntrials = numel(trialssort);
    
    xtemp = [0,1,1,0] -1 + xi*1.5;
    % ytemp = [0,0,1,1] + (i-1) * nt;
    for j = 0:ntrials-1
        ttr = ttr + 1;
        y = [ys(ttr), ys(ttr), ys(ttr+1), ys(ttr+1)];
        x = xtemp;
        % y = ytemp + j;
        patch('XData',x,'YData',y,'FaceColor',colors(colids(j+1),:))
        hold on
        mx = mean(x(1:2));
        my = mean(y(2:3));
        % text(mx,my,chans{i+1},'HorizontalAlignment','center','FontSize',6)
        text(mx,my,sprintf('%d',trialssort(j+1)),'HorizontalAlignment','center','FontSize',8,'FontUnits','normalized')
        hold on
    end
    
    nt = numel(trialssort);
    nr = numel(trialsr);
%     rep = sprintf('%d / %d (%d%%) %s',nr,nt,round(nr/nt*100),trialnames{i});
    
    ntvec(i) = nt;
    nrvec(i) = nr;
    
end

% %%
% % colors = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880];
% % colors = lines()
% 
% compst = 1:size(EEG.icaweights,2);
% compsr = EEG.preproc.comps;
% compsk = compst(~ismember(compst,compsr));
% 
% colids = [1*ones(1,numel(compsr)), 2*ones(1,numel(compsk))];
% compssort = [sort(compsr), sort(compsk)];
% 
% n = numel(compssort);
% 
% xtemp = [1,2,2,1] + xi;
% ytemp = [0,0,1,1];
% for i = 0:n-1
%     x = xtemp;
%     y = ytemp + i;
%     patch('XData',x,'YData',y,'FaceColor',colors(colids(i+1),:))
%     hold on
%     mx = mean(x(1:2));
%     my = mean(y(2:3));
%     % text(mx,my,chans{i+1},'HorizontalAlignment','center','FontSize',6)
%     text(mx,my,sprintf('%d',compssort(i+1)),'HorizontalAlignment','center','FontSize',8,'FontUnits','normalized')
%     hold on
% end
% % xlim([0,3])
% % ylim([0,n])
% na = numel(compssort);
% ni = numel(compsr);

% % Option 1:
% % set(gca,'xtick',3,'xticklabel',{sprintf('Rej. comps%s%s','\newline',rep)})
% 

% rep = sprintf('%d / %d (%d%%)',nr,nt,round(nr/nt*100));
% trialnames = {'L','R'};

for i = 1:numel(trialnames)
trialnames{i} = regexprep(trialnames{i}, '_', '');
end

rep = '';
for i = 1:nnames
    nt = ntvec(i);
    nr = nrvec(i);
    rep = [rep, sprintf('%d / %d (%d%%) %s',nr,nt,round(nr/nt*100),trialnames{i}), '\newline'];
end

if nnames > 1
    nt = sum(ntvec);
    nr = sum(nrvec);
    rep = [rep, sprintf('%d / %d (%d%%) %s',nr,nt,round(nr/nt*100)), '\newline'];
end

xti = mx;
xlab = {sprintf('Rej. trials%s%s','\newline',rep)};

end