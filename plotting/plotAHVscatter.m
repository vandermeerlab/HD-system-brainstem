function plotAHVscatter(cfg_in, sd)

cfg_def.smooth = 1;
cfg = ProcessConfig2(cfg_def, cfg_in);

%% plot scatterplot
AHV_dt = median(diff(sd.AHV.tvec));

cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = AHV_dt;
cfg_Q.tvec_edges = sd.AHV.tvec(1): AHV_dt: sd.AHV.tvec(end);

cfg_in.doPlot = 0;
[tc_out] = getAHV_TC(cfg_in, sd);

figure
for iCell = 1:length(sd.S.t)
    Stouse.t{1} = sd.S.t{iCell};
    F = MakeQfromS(cfg_Q, Stouse);  % convert to FR
    F.data = F.data ./ cfg_Q.dt;
    F_idx = nearest_idx3(sd.AHV.tvec, F.tvec);
    AHV_F = F.data(:,F_idx);
    ymax = max(AHV_F);
    set(gca, 'TickDir', 'out', 'FontSize', 16)
    plot(sd.AHV.data, AHV_F, '.', 'MarkerSize', .5, 'Color', [0.8500 0.3250 0.0980]); hold on
    xlabel('AHV deg./sec')
    ylabel('Firing Rate (Hz)')
    set(gca, 'Ylim', [0 ymax], 'FontSize', 16)
    
    if cfg.smooth == 1
        plot(tc_out.binCenters, smoothdata(tc_out.tc(iCell,:)), 'LineWidth', 5, 'Color', 'k');
    else
        plot(tc_out.binCenters, tc_out.tc(iCell,:), 'LineWidth', 5, 'Color', 'k');
    end
    title(sd.fn{iCell})
    disp('press any key to continue to the next cell')
    pause
    clf
end
