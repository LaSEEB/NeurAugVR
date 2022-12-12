function c = erp_component(lf, erp, nsources, spacing)
%% Noise
% Get sources
erp_source = lf_get_source_spaced(lf, nsources, spacing);
% ersp_source = lf_get_source_nearest(lf, pos);

c = struct();
c.signal = {erp};
c.source = erp_source;

% c_noise = utl_create_component(noise_source, noise, lf);
c = utl_check_component(c, lf);
end