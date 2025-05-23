function [tc_out] = getHD_TC(cfg_in, sd)
% JJS. 2023-08-24.
% Calculate and plot the Head Direction Tuning Curve
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

[csc_tsd, orientation, starttimeUNIX, endtimeUNIX, samplingrate, dt] = GetOrientationValues([]);

cfg_def = [];
cfg_def.doPlot = 1;
cfg_def.smooth = 0;
cfg_def.nBins = 60;
cfg_def.binEdges = {linspace(-180, 180, 101)};
cfg_def.occ_dt = dt;
cfg_def.minOcc = 200;  % remember that Occ is measured in samples, not in seconds. Usually 5ms per sample, b/c the platform encoder sampling rate is 200Hz,

cfg_tc = ProcessConfig2(cfg_def, cfg_in);


tc_out = TuningCurves(cfg_tc, sd.S, orientation);

if cfg_tc.doPlot
    for iCell = 1:length(sd.S.t)
        % Add Tuning Curve
        if cfg_tc.smooth
            plot(tc_out.binCenters, smoothdata(tc_out.tc(iCell,:)), 'LineWidth', 3, 'Color', 'k'); 
        else
            plot(tc_out.binCenters, tc_out.tc(iCell,:), 'LineWidth', 3, 'Color', 'k');
        end
        title(sd.S.label{iCell}); 
        xlabel('Head Direction (degrees)')
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