function [AHV_tsd] = Get_AHV(cfg_in)
%2020-03-12. JJS. Pulls out the encoder data and calculates AHV. 

cfg_def.subsample_factor = 10;

cfg = ProcessConfig(cfg_def, cfg_in);

[csc_tsd, hd_tsd, samplingrate, dt] = GetOrientationValues([]); %#ok<ASGLU>
orientationtouse = downsampleOrientationValues(hd_tsd, cfg.subsample_factor);

window = 0.1;
postsmoothing = .05;

tic 
AHV = dxdt(orientationtouse.tvec, orientationtouse.data, 'window', window, 'postsmoothing', postsmoothing);
toc

AHV = -AHV; % THIS STEP IS NECESSARY BECAUSE dxdt GIVES VALUES THAT ARE CORRECT, BUT WITH A SIGN FLIP.

AHV_tsd = tsd(orientationtouse.tvec, AHV);

end



