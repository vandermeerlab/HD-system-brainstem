% plot_EV_TC_from_workspace.m  (v2.0)
%
% Plots EV tuning curve with raw data scatter for a single neuron.
% Uses compute_slow_phase_EV for the EV pipeline (Order B).
%
% No longer requires data_outC or getSlowPhaseData2.
%
% Style matches getSlowPhaseTC3 / plot_AHVtuningMegaplot7:
%   - Grey dots [.8 .8 .8], MarkerSize 0.5
%   - TuningCurves (mvdm), smoothdata on TC values for well-occupied bins
%
% REQUIRES in workspace:
%   sd          - session data struct (from LoadSessionData([]))
%   iCell       - index of cell to plot (integer)
%
% USAGE:
%   sd = LoadSessionData([]);
%   iCell = 1;
%   plot_EV_TC_from_workspace
%
%                   - Fixed nearest_idx3 arg order
%                   - Interpolation within contiguous segments
% =========================================================================

%% ---- Settings ------------------------------------------------------------
FontSize    = 28;
MarkerSize  = 4;
LineWidth   = 5;
greyColor   = [0.45 0.45 0.45];
occthresh   = 0.5;                   % seconds of occupancy to include TC bin
binEdges    = linspace(-40, 40, 81); % deg/s
minOcc      = 0.5;                   % seconds (passed to TuningCurves)
gausswin_sd = 0.05;                  % s, Gaussian FR smoothing
%% -------------------------------------------------------------------------

SSN = HD_GetSSN;
myCell = SelectTS([], sd.S, iCell);

%% Compute slow-phase EV via standalone function
[ev_out, ev_params] = compute_slow_phase_EV();

%% Unpack interpolated EV
tvec_cont   = ev_out.ev_interp.tvec;
ev_cont     = ev_out.ev_interp.data;
ev_dt       = ev_out.ev_dt;

% Ensure column tvec, row data
if size(tvec_cont, 1) == 1; tvec_cont = tvec_cont'; end
if size(ev_cont, 1) ~= 1;  ev_cont = ev_cont'; end

%% Compute firing rate on continuous timebase
cfg_Q               = [];
cfg_Q.smooth        = 'gauss';
cfg_Q.gausswin_sd   = gausswin_sd;
cfg_Q.dt            = ev_dt;
cfg_Q.tvec_edges    = tvec_cont(1) : ev_dt : tvec_cont(end) + ev_dt;
F_ev                = MakeQfromS(cfg_Q, myCell);
F_ev.data           = F_ev.data ./ cfg_Q.dt;

% Align FR to continuous tvec
% nearest_idx3(A, B) → for each element of A, find nearest index into B
F_idx_ev    = nearest_idx3(tvec_cont, F_ev.tvec);
F_idx_ev    = min(F_idx_ev, size(F_ev.data, 2));  % guard against edge overrun
EV_FR_cont  = F_ev.data(1, F_idx_ev);

% Scatter: non-NaN slow-phase samples only
valid_mask  = ~isnan(ev_cont);
ev_scatter  = ev_cont(valid_mask);
fr_scatter  = EV_FR_cont(valid_mask);

%% Build tsd for TuningCurves (needs column tvec, row data on uniform grid)
% ev_cont is already 1×N row after orientation check above
ev_tsd = tsd(tvec_cont, ev_cont);

%% Compute EV tuning curve
cfg_tc          = [];
cfg_tc.nBins    = length(binEdges) - 1;
cfg_tc.binEdges = {binEdges};
cfg_tc.minOcc   = minOcc;
cfg_tc.occ_dt   = ev_dt;
tc_EV           = TuningCurves(cfg_tc, myCell, ev_tsd);

good_bins   = tc_EV.occ_hist > occthresh;
tc_smooth   = smoothdata(tc_EV.tc(1, good_bins));

%% Figure 1: Diagnostic — interpolated EV time series (first 60s)
fig_diag = figure('Position', [50 50 1600 400]);
plot(tvec_cont, ev_cont, 'k', 'LineWidth', 1);
xlabel('Time (s)',    'FontSize', FontSize - 4);
ylabel('EV (deg/s)', 'FontSize', FontSize - 4);
title(sprintf('%s — Slow-phase EV (Order B, first 60s)', SSN), ...
    'FontSize', FontSize - 4, 'Interpreter', 'none');
xlim([tvec_cont(1)  min(tvec_cont(1)+60, tvec_cont(end))]);
set(gca, 'FontSize', FontSize - 8, 'TickDir', 'out');
fig_diag.UserData.created_by   = 'plot_EV_TC_from_workspace';
fig_diag.UserData.version      = '2.0';
fig_diag.UserData.date_created = datestr(now, 'yyyy-mm-dd');
fig_diag.UserData.ev_params    = ev_params;

%% Figure 2: EV tuning curve
fig_tc = figure('Position', [50 500 900 700]); hold on;
set(gca, 'TickDir', 'out', 'FontSize', FontSize);

plot(ev_scatter, fr_scatter, '.', 'MarkerSize', MarkerSize, 'Color', greyColor);
plot(tc_EV.binCenters(good_bins), tc_smooth, 'k', 'LineWidth', LineWidth);

ymax_ev = max(tc_smooth) * 1.2;
if isempty(ymax_ev) || ~isfinite(ymax_ev) || ymax_ev <= 0
    ymax_ev = max(fr_scatter) * 1.2;
end
ylim([0 ymax_ev]);
xlim([binEdges(1) binEdges(end)]);
xlabel('Eye velocity (deg/s)', 'FontSize', FontSize);
ylabel('FR (Hz)',               'FontSize', FontSize);
title(sprintf('%s  cell %d — EV TC (Order B)', SSN, iCell), ...
    'FontSize', FontSize, 'Interpreter', 'none');

c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
text(0.15, 0.92, 'CW',  'Units', 'normalized', 'FontSize', FontSize - 2);
text(0.70, 0.92, 'CCW', 'Units', 'normalized', 'FontSize', FontSize - 2);

fig_tc.UserData.created_by   = 'plot_EV_TC_from_workspace';
fig_tc.UserData.version      = '2.0';
fig_tc.UserData.date_created = datestr(now, 'yyyy-mm-dd');
fig_tc.UserData.session      = SSN;
fig_tc.UserData.iCell        = iCell;
fig_tc.UserData.ev_params    = ev_params;
fig_tc.UserData.occthresh_s  = occthresh;

fprintf('EV TC plotted: %d good bins out of %d\n', sum(good_bins), length(good_bins));
fprintf('Valid scatter points: %d\n', sum(valid_mask));

%% Version History
% 2026-03-17  v2.0
%                   - Replaced data_outC dependency with compute_slow_phase_EV
%                   - Eliminates getSlowPhaseData2 dependency entirely
%                   - Order B pipeline now includes stationary-period removal
%                   - occthresh now in seconds (matching TuningCurves occ_hist
%                     units), was previously comparing seconds vs sample count
%                   - Added fallback for ymax_ev if tc_smooth is empty
% 2026-03-17  v1.1
