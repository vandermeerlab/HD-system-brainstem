function []  = saccadePETH(sd, cfg_in) 
% JJS. 2023-04-27. 
% This function pulls out the spike times in a window around temporal and nasal saccades

cfg_def.LineWidth = 3;
cfg_def.doPlot = 0;
cfg_def.window = [-.2 .2];
cfg_def.dt = 0.01;
cfg_out = ProcessConfig2(cfg_def, cfg_in);

[outputS, ~, ~, outputIT, ~] = SpikePETH_either(cfg_out, sd.S, sd.nasalSaccades);
[mn3, edges] = histcounts(outputS, outputIT);
plot(edges(1:end-1), mn3/cfg_out.dt/length(sd.nasalSaccades), 'LineWidth', cfg_out.LineWidth);
amax = max( mn3/cfg_1.dt/length(sd.nasalSaccades));
set(gca, 'FontSize', FontSize)

hold on
[outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, S, temporal_timestamps_MOVING);
[mt3, edges] = histcounts(outputS_t, outputIT_t);
plot(edges(1:end-1), mt3/cfg_1.dt/length(temporal_timestamps_MOVING), 'LineWidth', LineWidth);
bmax = max(mt3/cfg_1.dt/length(temporal_timestamps_MOVING));
allmax = max([amax bmax]);
title('Sacc. peth MOVING')
ylabel('FR (Hz)', 'FontSize', FontSize)
set(gca, 'FontSize', FontSize)
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
axis([c(1) c(2) c(3) allmax])