function [AHV_tsd, tc_out] = AHV_tuning(S, cfg_in)
%2020-03-12. JJS. Pulls out the encoder data and calculates AHV. Creates an tuning curve of AHV from velocity and spikes.

cfg_def.subsample_factor = 10;

cfg = ProcessConfig(cfg_def, cfg_in);

[csc_tsd, hd_tsd, samplingrate, dt] = GetOrientationValues([]); %#ok<ASGLU>
hd_ss_tsd = downsampleOrientationValues(hd_tsd, cfg.subsample_factor);
AHV_tsd = GetAHV_values(hd_ss_tsd);

% cfg = [];
% cfg.uint = '64';
% S = LoadSpikes(cfg);

cfg_tc = [];
cfg_tc.nBins = 100;
cfg_tc.binEdges = {linspace(-200, 200, 101)};
cfg_tc.occ_dt = median(diff(AHV_tsd.tvec));
cfg_tc.minOcc = 10;
tc_out = TuningCurves(cfg_tc, S, AHV_tsd);

end
