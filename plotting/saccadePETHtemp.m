function []  = saccadePETHtemp(sd, cfg_in)
% JJS. 2023-04-27.
% This function pulls out the spike times in a window around temporal and nasal saccades


cfg_def.smooth = 0;
cfg_def.FontSize = 16;
cfg_def.LineWidth = 3;
cfg_def.doPlot = 0;
cfg_def.window = [-.5 .5];
cfg_def.dt = 0.01;
cfg_def.doPlot = 1;
cfg_out = ProcessConfig2(cfg_def, cfg_in);

for iCell = 1:length(sd.S.t)
    Stouse.t{1} = sd.S.t{iCell};
    [outputS, ~, ~, outputIT, ~] = SpikePETHvdm(cfg_out, Stouse, sd.nasalSaccades); hold on
    subplot(2,1,1)
    set(gca, 'FontSize', cfg_out.FontSize)
    subplot(2,1,2)
    [m, edges] = histcounts(outputS, outputIT);
    if cfg_out.smooth
        FRnasal = smoothdata(m/cfg_out.dt/length(sd.nasalSaccades));
    else
        FRnasal = m/cfg_out.dt/length(sd.nasalSaccades);
    end
    
    
    if cfg_out.doPlot       
        plot(edges(1:end-1), FRnasal, 'LineWidth', cfg_out.LineWidth);
        amax = FRnasal;
        set(gca, 'FontSize', cfg_out.FontSize)
        disp('press any key to continue')
        pause
        clf
    end
end
% hold on
% [outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, S, temporal_timestamps_MOVING);
% [mt3, edges] = histcounts(outputS_t, outputIT_t);
% plot(edges(1:end-1), mt3/cfg_1.dt/length(temporal_timestamps_MOVING), 'LineWidth', LineWidth);
% bmax = max(mt3/cfg_1.dt/length(temporal_timestamps_MOVING));
% allmax = max([amax bmax]);
% title('Sacc. peth MOVING')
% ylabel('FR (Hz)', 'FontSize', FontSize)
% set(gca, 'FontSize', FontSize)
% c = axis;
% line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
% axis([c(1) c(2) c(3) allmax])