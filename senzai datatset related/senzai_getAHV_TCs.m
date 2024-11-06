function [tc_out_AHV] = senzai_getAHV_TCs(HD_neuronsToUse, S, AHVtsd, varargin)
%2024-10-31. JJS. This function calculates and plots the AHV tuning curves for all neurons
%   Inputs:
%               HD_neuronsToUse     - list of HD neurons for this mouse
%               S                   - structure with S.t{nNeurons} spike trains
%               AHVtsd              - tsd of angular head velocity


doPlot = 1;
doSmooth = 0;
process_varargin(varargin)
%% Calculate AHV tuning curves
AHV_dt = median(diff(AHVtsd.tvec));
cfg_tc = [];
cfg_tc.nBins = 67; % 67 bins from -200 to 200 works out to roughly 6 degree/s bins.
cfg_tc.binEdges = {linspace(-200, 200, cfg_tc.nBins)};
cfg_tc.occ_dt = median(diff(AHVtsd.tvec));
cfg_tc.minOcc = 100;  % remember that Occ is measured in samples (usually 20ms per sample), not in seconds
tc_out_AHV = TuningCurves(cfg_tc, S, AHVtsd);

if doPlot
    for iNeuron = 1:length(S.t)
        clf;
        if ismember(iNeuron, HD_neuronsToUse)
            colortouse = 'r';
        else
            colortouse = 'b';
        end
        if doSmooth
            plot(tc_out_AHV.binCenters, smoothdata(tc_out_AHV.tc(iNeuron,:)), 'Color', colortouse, 'LineWidth', 8)
        else
            plot(tc_out_AHV.binCenters, tc_out_AHV.tc(iNeuron,:), '.', 'Color', colortouse, 'MarkerSize', 20)
        end
        title(num2str(iNeuron))
        c = axis;
        axis([c(1) c(2) 0 c(4)]);
        xlabel('AHV (deg/s)')
        ylabel('FR (Hz)')
        set(gca, 'FontSize', 26)
        pause
    end
end

% if doPlot
%     %% #1 plot scatterplot
%     p = subtightplot(6,6,1, [tightX tightY]);
%     cfg_Q = [];
%     cfg_Q.smooth = 'gauss';
%     cfg_Q.gausswin_sd = 0.05;
%     cfg_Q.dt = AHV_dt;
%     cfg_Q.tvec_edges = AHVtsd.tvec(1):AHV_dt:AHVtsd.tvec(end);
%     F = MakeQfromS(cfg_Q, S); % convert to FR
%     F.data = F.data ./ cfg_Q.dt;
%     F_idx = nearest_idx3(AHVtsd.tvec, F.tvec);
%     AHV_F = F.data(:,F_idx);
%     ymax = max(AHV_F);
%     set(gca, 'TickDir', 'out', 'FontSize', FontSize)
%     plot(AHVtsd.data, AHV_F, '.', 'MarkerSize', .5); hold on
%     set(gca, 'Ylim', [0 ymax], 'FontSize', FontSize)
%     
%     % Add Tuning Curve
%     plot(tc_out_AHV.binCenters, tc_out_AHV.tc, 'LineWidth', LineWidth, 'Color', 'k');
%     ylabel('FR (Hz)', 'FontSize', FontSize)
%     set(groot, 'DefaultLegendInterpreter', 'none')
%     title('AHV Tuning Curve')
%     h = get(gca, 'XLim');
%     % text(.75*h(1), 10, 'CW', 'FontSize', 12)
%     % text(.5*h(2), 10, 'CCW', 'FontSize', 12)
%     text(NaN, NaN, 'CW', 'FontSize', insetText, 'Units', 'normalized', 'Position', [.25 .85 0])
%     text(NaN, NaN, 'CCW', 'FontSize', insetText, 'Units', 'normalized', 'Position', [.6 .85 0])
%     p.XAxisLocation = 'top';
%     c = axis;
%     line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
% end


end

