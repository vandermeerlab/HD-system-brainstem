function [] = plotAHVtuningplussaccadePETHs(varargin)
% 4-2021. JJS.
% Plot FR vs. AHV (scatterplot) with saccade PETHs overlaid in thick lines, for comparison.
FontSize = 13;
process_varargin(varargin);

[S] = LoadSpikesJeff;
% get AHV Tuning Curve
cfg_AHV = [];
cfg_AHV.subsample_factor = 10;
[AHV_tsd, tc_out] = AHV_tuning(S, cfg_AHV);
AHV_dt = median(diff(AHV_tsd.tvec));

[FRxBinT, FRxBinN, FRxBinTnorm, FRxBinNnorm, TnormSmooth, NnormSmooth, outputIT, cfg] = makeSaccadeHeatPlot;
for iCell = 1:length(S.t)
    %%  plot scatterplot
    figure(iCell)
    yyaxis left
    cfg_Q = [];
    cfg_Q.smooth = 'gauss';
    cfg_Q.gausswin_sd = 0.05;
    cfg_Q.dt = AHV_dt;
    F = MakeQfromS(cfg_Q, S); % convert to FR
    % convert to FR
    F.data = F.data ./ cfg_Q.dt;
    % find FR corresponding to each AHV sample
    F_idx = nearest_idx3(AHV_tsd.tvec, F.tvec);
    AHV_F = F.data(:,F_idx);
    ymax = max(AHV_F(iCell,:));
    
    set(gca, 'TickDir', 'out', 'FontSize', FontSize)
    plot(AHV_tsd.data, AHV_F(iCell,:), '.', 'MarkerSize', .5);
    set(gca, 'Ylim', [0 ymax], 'FontSize', FontSize)
    xlabel('AHV (deg./sec)', 'FontSize', FontSize)
    ylabel('FR (Hz)', 'FontSize', FontSize)
    title('Scatterplot')
    h = get(gca, 'XLim');
    text(.75*h(1), 10, 'CW', 'FontSize', 12)
    text(.5*h(2), 10, 'CCW', 'FontSize', 12)
    hold on
    
    yyaxis right
    plot(FRxBinT)
end
