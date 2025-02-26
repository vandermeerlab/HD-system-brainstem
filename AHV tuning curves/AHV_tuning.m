function [AHV_tsd, tc_out, cfg_TC] = AHV_tuning(cfg_AHV, S)
%2020-03-12. JJS. Pulls out the encoder data and calculates AHV. Creates an tuning curve of AHV from velocity and spikes.

cfg_def = [];
cfg_def.subsample_factor = 10;
cfg_def.nBins = 67; % 67 bins from -200 to 200 works out to roughly 6 degree/s bins. 
cfg_def.binEdges = {linspace(-200, 200, cfg_def.nBins)};
cfg_def.minOcc = 100;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds

cfg_to_use = ProcessConfig(cfg_def, cfg_AHV);

[csc_tsd, hd_tsd, samplingrate, dt] = GetOrientationValues([]); %#ok<ASGLU>
hd_ss_tsd = downsampleOrientationValues(hd_tsd, cfg_to_use.subsample_factor);
AHV_tsd = GetAHV_values(hd_ss_tsd);

% cfg = [];
% cfg.uint = '64';
% S = LoadSpikes(cfg);

% MvdM: this function will be more modular if you don't do the below here.
% I can see a number of use cases where you want to get AHV but don't care
% about tuning curves.

cfg_TC.occ_dt = median(diff(AHV_tsd.tvec));
tc_out = TuningCurves(cfg_to_use, S, AHV_tsd);

end
