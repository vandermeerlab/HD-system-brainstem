function [tc_out] = getAHV_TC(cfg_in, sd)
% JJS. 2023-04-27.
% Calculate and plot the AHV tuning curve
% input:   sd - session data structure with spike trains S and tsd of angular head velocity
% output:  tc_out - a structure
%               .usr.tc - nCells x mBins
%               .usr.occ_hist  -  occupancy in each bin
%               .usr.spk_hist  -  number of spikes in each AHV bin
%               .usr.good_idx  - bins with at least the minOcc # of samples
%               .usr.no_occ_idx  - bins with NaN values
%               .usr.binEdges
%               .usr.binCenters
%               .usr.nBins  - number of bins

% HELP
% Set 'doPlot'=1 to plot each TC. The default setting will smooth the data. 
% One parameter that you may want to change is minOCC. This value determines the minimum number of samples (each 5 ms) needed in a given AHV bin to include that data. 
% minOcc of 100 requires only 0.5 s of occupancy to include those values. I usually set this to 200 to include at least 1 s of data. But one may want to set it higher.

cfg_def = [];
cfg_def.doPlot = 1;
cfg_def.smooth = 1;
cfg_def.nBins = 100;
cfg_def.binEdges = {linspace(-200, 200, 101)};
cfg_def.occ_dt = median(diff(sd.AHV.tvec));
cfg_def.minOcc = 200;  % remember that Occ is measured in samples, not in seconds. Usually 5ms per sample, b/c the platform encoder sampling rate is 200Hz,

cfg_tc = ProcessConfig2(cfg_def, cfg_in);

tc_out = TuningCurves(cfg_tc, sd.S, sd.AHV);

if cfg_tc.doPlot
    for iCell = 1:length(sd.S.t)
        % Add Tuning Curve
        if cfg_tc.smooth
            plot(tc_out.binCenters, smoothdata(tc_out.tc(iCell,:)), 'LineWidth', 3, 'Color', 'k');
        else
            plot(tc_out.binCenters, tc_out.tc(iCell,:), 'LineWidth', 3, 'Color', 'k');
        end
        xlabel('AHV (degrees/sec)')
        ylabel('FR (Hz)')
        set(groot, 'DefaultLegendInterpreter', 'none')
        c = axis;
        line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
        set(gca, 'FontSize', 16)
        disp('press any key')
        pause
        clf
    end
end
