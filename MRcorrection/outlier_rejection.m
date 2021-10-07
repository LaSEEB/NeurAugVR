function [newdata, mavec, mivec] = outlier_rejection(data, method)

    newdata = zeros(size(data));
    mavec = zeros(size(data,1), 1);
    mivec = zeros(size(data,1), 1);
    switch method
        case 'capmean'
            for c = 1:size(data,1)
                chan = c;
                signal = data(chan, :);
                msig = mean(signal);
                stdsig = std(signal);
                ma = msig + 4*stdsig;
                mi = msig - 4*stdsig;

                out_ind = find((signal > ma) | (signal < mi));
                newsig = signal;
                win = 1;
                for k = 1:length(out_ind)
                    ind = out_ind(k);
                    if (ind > 1) & (ind < length(signal)-1)
                        sig_win = newsig(ind-win:ind+win);
                        sig_win(sig_win > ma) = ma;
                        sig_win(sig_win < mi) = mi;
                        val_ind = mean(sig_win);
                        newsig(ind) = val_ind;
                    elseif newsig(ind) > ma
                        newsig(ind) = ma;
                    elseif newsig(ind) < mi
                        newsig(ind) = mi;
                   end
                end
                newdata(c, :) = newsig;
                mavec(c,:) = ma;
                mivec(c,:) = mi;
            end
            
            
        case 'mean'
            
            for c = 1:size(data,1)
                chan = c;
                signal = data(chan, :);
                msig = mean(signal);
                stdsig = std(signal);
                ma = msig + 4*stdsig;
                mi = msig - 4*stdsig;

                out_ind = find((signal > ma) | (signal < mi));

                newsig = signal;
                win = 2;
                for k = 1:length(out_ind)
                    ind = out_ind(k);
                    if (ind > win) & (ind < length(signal)-win)
                        ind = out_ind(k);
                        sig_win = newsig(ind-win:ind+win);
                        val_ind = mean(sig_win);
                        newsig(ind) = val_ind;
                    elseif newsig(ind) > ma
                        newsig(ind) = ma;
                    elseif newsig(ind) < mi
                        newsig(ind) = mi;
                   end
                 end
                 newdata(c, :) = newsig;
                 mavec(c,:) = ma;
                 mivec(c,:) = mi;
            end
        otherwise
            contine
    end
end
