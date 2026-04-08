function diagnose_desaccade_order(varargin)
% diagnose_desaccade_order  Compare two EV filtering pipeline orders on a single session.
%
%   Order A (matches getSlowPhaseData2): desaccade -> medfilt1 -> smoothdata -> threshold
%   Order B (Jeff proposed):             desaccade -> threshold -> medfilt1
%
%   Plots a time-series snippet of the resulting slow-phase EV traces so
%   you can visually compare which order produces cleaner slow-phase data.
%
%   USAGE:
%       diagnose_desaccade_order()                        % uses defaults
%       diagnose_desaccade_order('sess_dir', 'C:\...')    % override session
%       diagnose_desaccade_order('snippet_start', 120)    % set window start
%
%   OUTPUTS:
%       Saves .png to ___analysis_output\diagnose_desaccade_<timestamp>\
%

%% Defaults
sess_dir        = 'C:\Jeff\U01\datatouse\M222\M222-2022-01-21';
snippet_start   = [];    % leave empty to auto-pick (5s into de-saccaded data)
snippet_dur     = 30;    % seconds to display
saccade_pre     = 0.1;
saccade_post    = 0.1;
medfilt_win     = 11;    % samples, matches getSlowPhaseData2
vel_thresh      = 35;    % deg/s, matches getSlowPhaseData2
FontSize        = 28;
save_root       = 'C:\Jeff\U01\datatouse\___analysis_output';

%% Process varargin
process_varargin(varargin);

%% Store params
params = struct();
params.function_name  = 'diagnose_desaccade_order';
params.date_run       = datestr(now, 'yyyy-mm-dd');
params.sess_dir       = sess_dir;
params.snippet_dur    = snippet_dur;
params.saccade_pre    = saccade_pre;
params.saccade_post   = saccade_post;
params.medfilt_win    = medfilt_win;
params.vel_thresh     = vel_thresh;

%% Load session data
pushdir(sess_dir);
SSN = HD_GetSSN;
fprintf('Session: %s\n', SSN);

sacc_data = load(strcat(SSN, '-saccades-edited.mat'), ...
    'combinedSaccades', 'diffHdeg');
diffHdeg_raw  = sacc_data.diffHdeg;
combinedSacc  = sacc_data.combinedSaccades;

% Ensure correct orientation
if size(diffHdeg_raw.tvec, 1) == 1
    diffHdeg_raw.tvec = diffHdeg_raw.tvec';
end
if size(diffHdeg_raw.data, 1) ~= 1
    diffHdeg_raw.data = diffHdeg_raw.data';
end

%% Build saccade exclusion windows
sacc_times  = combinedSacc(~isnan(combinedSacc));
pre_sacc    = sacc_times - saccade_pre;
post_sacc   = sacc_times + saccade_post;

%% Step 1 (shared): de-saccade raw EV
[~, keep_mask] = restrict(diffHdeg_raw, pre_sacc, post_sacc);
is_desacc  = keep_mask == 0;
tvec_ds    = diffHdeg_raw.tvec(is_desacc);
ev_ds      = diffHdeg_raw.data(is_desacc);

%% Order A: desaccade -> medfilt1 -> smoothdata -> threshold
ev_A_filt         = medfilt1(ev_ds, medfilt_win);
ev_A_smooth       = smoothdata(ev_A_filt);
keep_A            = abs(ev_A_smooth) < vel_thresh;
tvec_A            = tvec_ds(keep_A);
ev_A_raw_kept     = ev_ds(keep_A);
ev_A_smooth_kept  = ev_A_smooth(keep_A);

%% Order B: desaccade -> threshold -> medfilt1
keep_B_thresh   = abs(ev_ds) < vel_thresh;
tvec_B          = tvec_ds(keep_B_thresh);
ev_B_raw_kept   = ev_ds(keep_B_thresh);
ev_B_filt_kept  = medfilt1(ev_B_raw_kept, medfilt_win);

%% Pick snippet window
if isempty(snippet_start)
    snippet_start = tvec_ds(1) + 5;
end
snippet_end = snippet_start + snippet_dur;

in_snip = @(tv) tv >= snippet_start & tv <= snippet_end;

mask_ds  = in_snip(tvec_ds);
mask_A   = in_snip(tvec_A);
mask_B   = in_snip(tvec_B);

fprintf('Snippet: %.1f - %.1f s\n', snippet_start, snippet_end);
fprintf('Order A: %d samples in snippet\n', sum(mask_A));
fprintf('Order B: %d samples in snippet\n', sum(mask_B));

%% Plot
fig_diag = figure('Position', [50 100 1800 900], 'Visible', 'on');
fig_diag.UserData.created_by   = 'diagnose_desaccade_order';
fig_diag.UserData.date_created = datestr(now, 'yyyy-mm-dd');
fig_diag.UserData.session      = SSN;
fig_diag.UserData.params       = params;

thresh_col = [0.8 0 0];
grey_col   = [0.5 0.5 0.5];

% --- Row 1: Raw de-saccaded EV (shared baseline) ---
ax1 = subplot(3, 1, 1); hold on;
plot(tvec_ds(mask_ds), ev_ds(mask_ds), 'Color', grey_col, 'LineWidth', 1.5);
yline( vel_thresh, '--', 'Color', thresh_col, 'LineWidth', 2);
yline(-vel_thresh, '--', 'Color', thresh_col, 'LineWidth', 2);
ylabel('EV (deg/s)', 'FontSize', FontSize);
title(sprintf('%s  —  Raw de-saccaded EV (shared input)', SSN), ...
    'FontSize', FontSize, 'Interpreter', 'none');
set(gca, 'FontSize', FontSize-4, 'TickDir', 'out');
xlim([snippet_start snippet_end]);

% --- Row 2: Order A ---
ax2 = subplot(3, 1, 2); hold on;
plot(tvec_A(mask_A), ev_A_raw_kept(mask_A), ...
    'Color', [0.6 0.6 1.0], 'LineWidth', 1.5);
plot(tvec_A(mask_A), ev_A_smooth_kept(mask_A), ...
    'b', 'LineWidth', 2.5);
yline( vel_thresh, '--', 'Color', thresh_col, 'LineWidth', 2);
yline(-vel_thresh, '--', 'Color', thresh_col, 'LineWidth', 2);
ylabel('EV (deg/s)', 'FontSize', FontSize);
title('Order A:  desaccade  \rightarrow  medfilt1  \rightarrow  smoothdata  \rightarrow  threshold   (getSlowPhaseData2)', ...
    'FontSize', FontSize-4);
legend({'raw slow-phase EV', 'smoothed (used for threshold)'}, ...
    'FontSize', FontSize-6, 'Location', 'northeast');
set(gca, 'FontSize', FontSize-4, 'TickDir', 'out');
xlim([snippet_start snippet_end]);

% --- Row 3: Order B ---
ax3 = subplot(3, 1, 3); hold on;
plot(tvec_B(mask_B), ev_B_raw_kept(mask_B), ...
    'Color', [0.4 0.8 0.4], 'LineWidth', 1.5);
plot(tvec_B(mask_B), ev_B_filt_kept(mask_B), ...
    'Color', [0.0 0.5 0.0], 'LineWidth', 2.5);
yline( vel_thresh, '--', 'Color', thresh_col, 'LineWidth', 2);
yline(-vel_thresh, '--', 'Color', thresh_col, 'LineWidth', 2);
ylabel('EV (deg/s)', 'FontSize', FontSize);
xlabel('Time (s)', 'FontSize', FontSize);
title('Order B:  desaccade  \rightarrow  threshold  \rightarrow  medfilt1   (Jeff proposed)', ...
    'FontSize', FontSize-4);
legend({'raw slow-phase EV', 'medfilt1 after threshold'}, ...
    'FontSize', FontSize-6, 'Location', 'northeast');
set(gca, 'FontSize', FontSize-4, 'TickDir', 'out');
xlim([snippet_start snippet_end]);

linkaxes([ax1 ax2 ax3], 'xy');

%% Save
timestamp_str = datestr(now, 'yyyy-mm-dd_hh-MM-AM');
out_folder = fullfile(save_root, ['diagnose_desaccade_' timestamp_str]);
mkdir(out_folder);

png_path = fullfile(out_folder, 'diagnose_desaccade_order.png');
print(fig_diag, png_path, '-dpng', '-r150');
save(fullfile(out_folder, 'diagnose_desaccade_params.mat'), 'params');
fprintf('Saved to: %s\n', out_folder);

popdir;

%% Version History
% 2026-03-17  v1.0 Initial version.
