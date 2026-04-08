% plot_EV_TC_from_workspace.m
%
% Plots EV tuning curve with raw data scatter from workspace variables.
% Style matches plot_AHVtuningMegaplot7 pupil TC panel exactly:
%   - Grey dots [.8 .8 .8], MarkerSize 0.5
%   - TC via TuningCurves, smoothdata applied only to bins where occ_hist > occthresh
%   - TickDir out
%
% REQUIRES in workspace:
%   sd          - session data struct (from LoadSessionData([]))
%   data_outC   - slow-phase data struct (from getSlowPhaseData2([]))
%                 uses data_outC.horiz_eye_vel (desaccaded + thresholded + filtered)
%   iCell       - index of cell to plot (integer)
%
% USAGE:
%   iCell = 1;
%   plot_EV_TC_from_workspace

%% ---- Settings (edit here) ------------------------------------------------
FontSize    = 28;
MarkerSize  = 0.5;
LineWidth   = 5;
occthresh   = 0.5;                     % min occupancy to include bin in TC
binEdges    = linspace(-40, 40, 81);   % deg/s
minOcc      = 25;                      % samples
gausswin_sd = 0.05;                    % s, for Q matrix smoothing
%% -------------------------------------------------------------------------

SSN = HD_GetSSN;

%% Select neuron
myCell = SelectTS([], sd.S, iCell);

%% Get slow-phase EV tsd from data_outC
ev_slow = data_outC.horiz_eye_vel;   % tsd: desaccaded + slow-phase cutoff

% Ensure correct orientation
if size(ev_slow.tvec, 1) == 1
    ev_slow.tvec = ev_slow.tvec';
end
if size(ev_slow.data, 1) ~= 1
    ev_slow.data = ev_slow.data';
end

%% Compute firing rate on slow-phase timebase
ev_dt = median(diff(ev_slow.tvec));

cfg_Q               = [];
cfg_Q.smooth        = 'gauss';
cfg_Q.gausswin_sd   = gausswin_sd;
cfg_Q.dt            = ev_dt;
cfg_Q.tvec_edges    = ev_slow.tvec(1):ev_dt:ev_slow.tvec(end);
F_ev                = MakeQfromS(cfg_Q, myCell);
F_ev.data           = F_ev.data ./ cfg_Q.dt;

F_idx_ev    = nearest_idx3(ev_slow.tvec, F_ev.tvec);
EV_FR       = F_ev.data(:, F_idx_ev);

%% Compute EV tuning curve
cfg_tc          = [];
cfg_tc.nBins    = length(binEdges) - 1;
cfg_tc.binEdges = {binEdges};
cfg_tc.minOcc   = minOcc;
cfg_tc.occ_dt   = ev_dt;
tc_EV           = TuningCurves(cfg_tc, myCell, ev_slow);

%% Plot
figure; hold on;
set(gca, 'TickDir', 'out', 'FontSize', FontSize);

% Raw data cloud — grey, small dots (matches megaplot7 pupil TC panel)
plot(ev_slow.data, EV_FR(1,:), '.', 'MarkerSize', MarkerSize, 'Color', [.8 .8 .8]);

% TC overlay — smoothdata only on well-occupied bins (matches megaplot7)
good_bins = tc_EV.occ_hist > occthresh;
tc_smooth = smoothdata(tc_EV.tc(1, good_bins));
plot(tc_EV.binCenters(good_bins), tc_smooth, 'k', 'LineWidth', LineWidth);

% Formatting
ymax_ev = max(EV_FR(isfinite(EV_FR)));
ylim([0 ymax_ev * 1.05]);
xlim([binEdges(1) binEdges(end)]);
xlabel('Eye velocity (deg/s)', 'FontSize', FontSize);
ylabel('FR (Hz)',               'FontSize', FontSize);
title(sprintf('%s  cell %d', SSN, iCell), 'FontSize', FontSize, 'Interpreter', 'none');

c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
text(0.15, 0.92, 'CW',  'Units', 'normalized', 'FontSize', FontSize - 2);
text(0.70, 0.92, 'CCW', 'Units', 'normalized', 'FontSize', FontSize - 2);

% Metadata
fig_h = gcf;
fig_h.UserData.created_by   = 'plot_EV_TC_from_workspace';
fig_h.UserData.date_created = datestr(now, 'yyyy-mm-dd');
fig_h.UserData.session      = SSN;
fig_h.UserData.iCell        = iCell;
