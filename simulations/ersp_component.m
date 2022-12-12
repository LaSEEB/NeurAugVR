function c = ersp_component(lf, ersp, pos)

% Get sources
ersp_source = lf_get_source_nearest(lf, pos);

% Set activation pattern - ERSP
% (on top)


% Set components
c = struct();
c.signal = {ersp};         % ERSP class, defined above
c.source = ersp_source;    % obtained from the lead field, as above
% c.orientation = [-lf.chanlocs(22).Y, lf.chanlocs(22).X, lf.chanlocs(22).Z];
% c.orientation = [-1, 0, .75];

c = utl_check_component(c, lf);
% c = utl_shift_latency(c, delay);

end