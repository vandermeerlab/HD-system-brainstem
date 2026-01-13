%% #15 Laser PETH
% p = subtightplot(6,6,15, [tightX tightY]); hold on
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Get the correct timestamps for the different session events
% [start_time, stop_time, laser_on, laser_off, bit0, bit4, arraysize] = SortBrainstemEventLabels;
[start_time, stop_time, laser_on, laser_off, arraysize, difflaser] = SortBrainstemEventLabels3;

dur = mode(laser_off(1:arraysize) - laser_on(1:arraysize));

 sd = LoadSessionData([]);  
iCell = 2;
myCell = SelectTS([], sd.S, iCell);

%%
clf
cfg_laser.window = [-1.1 2];
cfg_laser.dt = 0.01;
cfg_laser.binsize = cfg_laser.dt; % used for gaussian kernal.  select a small bin size for good time resolution
cfg_laser.doPlot = 1;
cfg_laser.doRaster = 1;
cfg_laser.doBar = 0;
[~, ~, ~, ~, ~] = SpikePETH_either(cfg_laser, myCell, laser_on); hold on
c = axis;
rectangle(Position = [0, 0, dur, c(4)], FaceColor=[0 1 1], EdgeColor=[0 1 1])    % change the 3rd entry in rectangle to the diff of laser_off and laser_on
[outputS, outputT, outputGau, outputIT, cfg] = SpikePETH_either(cfg_laser, myCell, laser_on);  % this is a hack to get the background color 'in back'. otherwise it occludes the spikes.
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'w')
line([1 1], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'w')
set(gca, 'YTick', [])
set(gca, 'FontSize', FontSize)
title('Laser PETH', 'FontSize', FontSize)
set(gca, 'XTick', [-1 0 1 2])
hold on 
% add line for Ave Firing Rate, esp. when rasters are so dense that it turns black
cfg_laser.dt = 0.05;
cfg_laser.doPlot = 0;
cfg_laser.doRaster = 0;
cfg_laser.doBar = 0;
[outputS, outputT, outputGau, outputIT, cfg] = SpikePETH_either(cfg_laser, myCell, laser_on);
m = histc(outputS, outputIT);
yyaxis right
plot(outputIT,m/cfg.dt/length(laser_on), 'LineWidth', 6);