function plotAHVscatter(sd, varargin)

smooth = 1;
process_varargin(varargin);

%% plot scatterplot
AHV_dt = median(diff(sd.AHV.tvec));

cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = AHV_dt;
cfg_Q.tvec_edges = sd.AHV.tvec(1): AHV_dt: sd.AHV.tvec(end);

for iCell = 1:length(sd.S.t)
    Stouse.t{1} = sd.S.t{iCell}; 
    F = MakeQfromS(cfg_Q, Stouse);  % convert to FR
    F.data = F.data ./ cfg_Q.dt;
    F_idx = nearest_idx3(sd.AHV.tvec, F.tvec);
    AHV_F = F.data(:,F_idx);
    ymax = max(AHV_F);
    set(gca, 'TickDir', 'out', 'FontSize', 16)
    plot(sd.AHV.data, AHV_F, '.', 'MarkerSize', .5); hold on
    set(gca, 'Ylim', [0 ymax], 'FontSize', 16)
    
    cfg_in.doPlot = 0;
    [tc_out] = getAHV_TC(sd, cfg_in);
    if smooth == 1
        plot(tc_out.usr.binCenters, smoothdata(tc_out.tc(iCell,:)), 'LineWidth', 5, 'Color', 'k');
    else
        plot(tc_out.usr.binCenters, tc_out.tc(iCell,:), 'LineWidth', 5, 'Color', 'k');
    end
    disp('press any key to continue to the next cell')
    pause
    clf
end
