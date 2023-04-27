function [tc_out] = getAHV_TC(sd)
% JJS. 2023-04-27.
% input:   sd - session data structure with spike trains S and AHVtsd
% output:  tc_out  


cfg_tc = [];
cfg_tc.nBins = 100;
cfg_tc.binEdges = {linspace(-200, 200, 101)};
cfg_tc.occ_dt = median(diff(sd.AHV.tvec));
cfg_tc.minOcc = 100;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
tc_out = TuningCurves(cfg_tc, sd.S, sd.AHV);

end