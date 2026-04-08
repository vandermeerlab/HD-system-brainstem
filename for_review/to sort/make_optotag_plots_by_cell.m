%% Make AHV plot w/ raw data and TC average
clf
iCell = 2;
sd = LoadSessionData([]);
plot_AHV_TC_w_raw_data(sd, iCell)
box off
set(gca, 'FontSize', 18)
set(gca, 'TickDir', 'out');
xlabel('AHV (deg/s)')
ylabel('Firing Rate (Hz)')
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Color', 'k');


%% Saccade peth plotted as lines
clf
myCell = SelectTS([], sd.S, iCell);
[spkTimes_N, spkTimes_T, eventindexN, eventindexT, FR_N, FR_T, centers, edges, cfg_out]  = saccadePETH_one_neuron([], sd.nasalSaccades, sd.temporalSaccades, myCell);
c = axis;
axis([c(1) c(2) 0 c(4)]);
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Color', 'k');
box off
set(gca, 'FontSize', 18)
set(gca, 'TickDir', 'out');
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
legend('preferred', 'non-preferred')

%%  Laser peth as bar plot
clf
[start_time, stop_time, laser_on, laser_off, arraysize, difflaser] = SortBrainstemEventLabels3;
cfg_in.window = [-1 2];
[outputS, outputT, outputGau, outputIT, cfg] = make_laser_PETH(cfg_in, myCell, laser_on);

%% Second laser peth option
clf
plot_laser_PETH(sd, iCell)
box off
set(gca, 'FontSize', 18)
set(gca, 'TickDir', 'out');
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
title('')
c = axis;
axis([c(1) c(2) 0 c(4)]);