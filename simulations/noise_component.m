function c = noise_component(lf, noise, nsources, spacing)
%% Noise
% Get sources
noise_source = lf_get_source_spaced(lf, nsources, spacing);

c = struct();
c.signal = {noise};
c.source = noise_source;

% c_noise = utl_create_component(noise_source, noise, lf);
c = utl_check_component(c, lf);

end