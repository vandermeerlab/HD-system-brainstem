%% plot_laser_eye_peth_diagnostic_ver_1_0.m  — v1.0
%
% Plots laser-aligned eye position and eye velocity PETHs for the 7
% optotagged sessions WITHOUT absolute value. Shows whether laser
% stimulation drives consistent directional eye movements.
%
% Two figures:
%   1) Signed eye position PETH (baseline-subtracted, no abs)
%   2) Eye velocity PETH (derivative of eye position)
%
% Each figure: 7 per-session subplots + 1 population average subplot.
%
% Requires: nothing in workspace. Loads all data from disk.
%
% Version History:
%   2026-04-07  v1.0  Initial version.

%% =====================================================================
%  CONFIG
%  =====================================================================
data_root  = 'C:\Jeff\U01\datatouse';
save_dir   = 'C:\Jeff\U01\datatouse\___analysis_output';

opto_sessions = {'M271-2021-08-30', 'M282-2022-02-10', 'M282-2022-02-15', ...
    'M391-2022-12-07-1', 'M391-2022-12-07-2', 'M402-2022-11-02-2', ...
    'M444-2023-08-18'};

align_window = [-1 2];
dt           = 0.005;
baseline_win = [-1 -0.2];
smooth_win   = 5;
FontSize     = 16;
TitleFS      = 14;
LineWidth    = 3;

%% Output folder
run_date_str = datestr(now, 'yyyy-mm-dd');
run_time_str = datestr(now, 'hh-MM-SS_AM');
out_folder = fullfile(save_dir, sprintf('laser_eye_peth_diag_%s_%s', run_date_str, run_time_str));
mkdir(out_folder);

%% Time grid
eye_grid = align_window(1) : dt : align_window(2);
eye_grid = eye_grid(1:end-1) + dt / 2;
nBins    = length(eye_grid);
bl_mask  = eye_grid >= baseline_win(1) & eye_grid <= baseline_win(2);

% Velocity grid (one fewer bin due to diff)
vel_grid = eye_grid(1:end-1) + dt/2;
nBins_vel = length(vel_grid);

nSessions = length(opto_sessions);

%% =====================================================================
%  LOAD AND COMPUTE
%  =====================================================================
pos_means = NaN(nSessions, nBins);       % signed position, per session mean
vel_means = NaN(nSessions, nBins_vel);   % velocity, per session mean
n_trials_per_sess = NaN(nSessions, 1);
laser_durs = NaN(nSessions, 1);
session_ok = false(nSessions, 1);

for iS = 1:nSessions
    SSN = opto_sessions{iS};
    mouse_id = SSN(1:4);
    sess_path = fullfile(data_root, mouse_id, SSN);
    fprintf('=== %d/%d: %s ===\n', iS, nSessions, SSN);

    if contains(SSN, 'M271'); laser_durs(iS) = 0.5;
    else; laser_durs(iS) = 1.0; end

    try
        pushdir(sess_path);

        % Laser times
        [~, ~, laser_on_times] = SortBrainstemEventLabels3;
        if isnan(laser_on_times); fprintf('  No laser\n'); popdir; continue; end

        % Eye data
        sacc_file = strcat(SSN, '-saccades-edited.mat');
        if ~exist(sacc_file, 'file'); fprintf('  No saccades file\n'); popdir; continue; end
        eye_vars = load(sacc_file, 'tsdH', 'k');
        if ~isfield(eye_vars, 'tsdH'); fprintf('  No tsdH\n'); popdir; continue; end

        ep_tvec = eye_vars.tsdH.tvec(:);
        ep_data = eye_vars.tsdH.data(:);
        if isfield(eye_vars, 'k'); cal_k = eye_vars.k; else; cal_k = 1; end
        ep_deg = ep_data * cal_k;

        % Align each trial
        nTrials = length(laser_on_times);
        trial_pos = NaN(nTrials, nBins);

        for iT = 1:nTrials
            t0 = laser_on_times(iT);
            search_mask = ep_tvec >= (t0 + align_window(1) - 0.1) & ...
                          ep_tvec <= (t0 + align_window(2) + 0.1);
            if sum(search_mask) < 5; continue; end

            rel_t = ep_tvec(search_mask) - t0;
            rel_d = ep_deg(search_mask);
            interp_trace = interp1(rel_t, rel_d, eye_grid, 'linear', NaN);

            % Baseline subtract (signed, no abs)
            bl_val = mean(interp_trace(bl_mask), 'omitnan');
            if isnan(bl_val); continue; end
            trial_pos(iT, :) = interp_trace - bl_val;
        end

        valid = ~all(isnan(trial_pos), 2);
        trial_pos_valid = trial_pos(valid, :);
        n_valid = sum(valid);

        if n_valid < 2; fprintf('  <2 valid trials\n'); popdir; continue; end

        % Mean signed position
        pos_means(iS, :) = mean(trial_pos_valid, 1, 'omitnan');

        % Eye velocity: diff of each trial, then average
        trial_vel = diff(trial_pos_valid, 1, 2) / dt;  % deg/s
        vel_means(iS, :) = mean(trial_vel, 1, 'omitnan');

        n_trials_per_sess(iS) = n_valid;
        session_ok(iS) = true;
        fprintf('  %d valid trials\n', n_valid);
        popdir;
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        try popdir; catch; end
    end
end

%% Smooth
for iS = 1:nSessions
    if session_ok(iS)
        pos_means(iS,:) = smoothdata(pos_means(iS,:), 'gaussian', smooth_win);
        vel_means(iS,:) = smoothdata(vel_means(iS,:), 'gaussian', smooth_win);
    end
end

%% Population averages
n_ok = sum(session_ok);
pop_pos_mean = mean(pos_means(session_ok, :), 1, 'omitnan');
pop_pos_sem  = std(pos_means(session_ok, :), 0, 1, 'omitnan') / sqrt(n_ok);
pop_vel_mean = mean(vel_means(session_ok, :), 1, 'omitnan');
pop_vel_sem  = std(vel_means(session_ok, :), 0, 1, 'omitnan') / sqrt(n_ok);

%% =====================================================================
%  FIGURE 1: SIGNED EYE POSITION PETH
%  =====================================================================
fig_pos = figure('Position', [30 50 2000 600]);
fig_pos.UserData.date_created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
fig_pos.UserData.generating_function = 'plot_laser_eye_peth_diagnostic_ver_1_0';

color_line = [0.60 0.15 0.60];
color_shade = [0.80 0.60 0.80];

for iS = 1:nSessions
    ax = subplot(1, nSessions + 1, iS); hold(ax, 'on');
    if session_ok(iS)
        ld = laser_durs(iS);
        % Laser shading
        yl_est = max(abs(pos_means(iS,:))) * 1.3;
        if yl_est == 0; yl_est = 1; end
        fill(ax, [0 ld ld 0], [-yl_est -yl_est yl_est yl_est], ...
            [0.4 0.7 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
        plot(ax, eye_grid, pos_means(iS,:), '-', 'Color', color_line, 'LineWidth', LineWidth);
        yline(ax, 0, '-', 'Color', [.6 .6 .6], 'LineWidth', 0.5);
        ylim(ax, [-yl_est yl_est]);
        title(ax, sprintf('%s\nn=%d', opto_sessions{iS}, n_trials_per_sess(iS)), ...
            'FontSize', TitleFS, 'Interpreter', 'none');
    else
        text(ax, 0.5, 0.5, 'No data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', [.6 .6 .6]);
        title(ax, opto_sessions{iS}, 'FontSize', TitleFS, 'Interpreter', 'none');
    end
    xline(ax, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1);
    xlim(ax, align_window); box(ax, 'off'); set(ax, 'FontSize', 13);
    if iS == 1; ylabel(ax, 'Eye position (deg)', 'FontSize', FontSize); end
    if iS == 4; xlabel(ax, 'Time from laser onset (s)', 'FontSize', FontSize); end
end

% Population
ax_pop = subplot(1, nSessions + 1, nSessions + 1); hold(ax_pop, 'on');
fill(ax_pop, [0 1 1 0], [-5 -5 5 5], [0.4 0.7 1.0], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.15);
fx = [eye_grid fliplr(eye_grid)];
fy = [pop_pos_mean + pop_pos_sem, fliplr(pop_pos_mean - pop_pos_sem)];
if ~any(isnan(fy))
    fill(ax_pop, fx, fy, color_shade, 'EdgeColor', 'none', 'FaceAlpha', 0.35);
end
plot(ax_pop, eye_grid, pop_pos_mean, '-', 'Color', color_line, 'LineWidth', LineWidth);
yline(ax_pop, 0, '-', 'Color', [.6 .6 .6], 'LineWidth', 0.5);
xline(ax_pop, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1);
xlim(ax_pop, align_window); box(ax_pop, 'off'); set(ax_pop, 'FontSize', 13);
title(ax_pop, sprintf('Population\n(N=%d)', n_ok), 'FontSize', TitleFS, 'FontWeight', 'bold');
sgtitle('Laser — Signed Eye Position (no abs)', 'FontSize', 20, 'FontWeight', 'bold');

print(fig_pos, fullfile(out_folder, 'laser_signed_eye_position.png'), '-dpng', '-r300');
savefig(fig_pos, fullfile(out_folder, 'laser_signed_eye_position.fig'));

%% =====================================================================
%  FIGURE 2: EYE VELOCITY PETH
%  =====================================================================
fig_vel = figure('Position', [30 100 2000 600]);
fig_vel.UserData.date_created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
fig_vel.UserData.generating_function = 'plot_laser_eye_peth_diagnostic_ver_1_0';

color_vel = [0.15 0.50 0.15];
color_vel_shade = [0.55 0.80 0.55];

for iS = 1:nSessions
    ax = subplot(1, nSessions + 1, iS); hold(ax, 'on');
    if session_ok(iS)
        ld = laser_durs(iS);
        yl_est = max(abs(vel_means(iS,:))) * 1.3;
        if yl_est == 0; yl_est = 5; end
        fill(ax, [0 ld ld 0], [-yl_est -yl_est yl_est yl_est], ...
            [0.4 0.7 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
        plot(ax, vel_grid, vel_means(iS,:), '-', 'Color', color_vel, 'LineWidth', LineWidth);
        yline(ax, 0, '-', 'Color', [.6 .6 .6], 'LineWidth', 0.5);
        ylim(ax, [-yl_est yl_est]);
        title(ax, sprintf('%s\nn=%d', opto_sessions{iS}, n_trials_per_sess(iS)), ...
            'FontSize', TitleFS, 'Interpreter', 'none');
    else
        text(ax, 0.5, 0.5, 'No data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', [.6 .6 .6]);
        title(ax, opto_sessions{iS}, 'FontSize', TitleFS, 'Interpreter', 'none');
    end
    xline(ax, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1);
    xlim(ax, align_window); box(ax, 'off'); set(ax, 'FontSize', 13);
    if iS == 1; ylabel(ax, 'Eye velocity (deg/s)', 'FontSize', FontSize); end
    if iS == 4; xlabel(ax, 'Time from laser onset (s)', 'FontSize', FontSize); end
end

% Population
ax_pop2 = subplot(1, nSessions + 1, nSessions + 1); hold(ax_pop2, 'on');
fill(ax_pop2, [0 1 1 0], [-20 -20 20 20], [0.4 0.7 1.0], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.15);
fx_v = [vel_grid fliplr(vel_grid)];
fy_v = [pop_vel_mean + pop_vel_sem, fliplr(pop_vel_mean - pop_vel_sem)];
if ~any(isnan(fy_v))
    fill(ax_pop2, fx_v, fy_v, color_vel_shade, 'EdgeColor', 'none', 'FaceAlpha', 0.35);
end
plot(ax_pop2, vel_grid, pop_vel_mean, '-', 'Color', color_vel, 'LineWidth', LineWidth);
yline(ax_pop2, 0, '-', 'Color', [.6 .6 .6], 'LineWidth', 0.5);
xline(ax_pop2, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1);
xlim(ax_pop2, align_window); box(ax_pop2, 'off'); set(ax_pop2, 'FontSize', 13);
title(ax_pop2, sprintf('Population\n(N=%d)', n_ok), 'FontSize', TitleFS, 'FontWeight', 'bold');
sgtitle('Laser — Eye Velocity (signed)', 'FontSize', 20, 'FontWeight', 'bold');

print(fig_vel, fullfile(out_folder, 'laser_eye_velocity.png'), '-dpng', '-r300');
savefig(fig_vel, fullfile(out_folder, 'laser_eye_velocity.fig'));

%% Save params
diag_params = struct();
diag_params.function_name = 'plot_laser_eye_peth_diagnostic_ver_1_0';
diag_params.version       = '1.0';
diag_params.date_run      = datestr(now, 'yyyy-mm-dd HH:MM:SS');
diag_params.align_window  = align_window;
diag_params.dt            = dt;
diag_params.baseline_win  = baseline_win;
diag_params.smooth_win    = smooth_win;
diag_params.n_sessions_ok = n_ok;

save(fullfile(out_folder, 'laser_eye_peth_diagnostic.mat'), ...
    'pos_means', 'vel_means', 'pop_pos_mean', 'pop_pos_sem', ...
    'pop_vel_mean', 'pop_vel_sem', 'eye_grid', 'vel_grid', ...
    'n_trials_per_sess', 'session_ok', 'opto_sessions', 'diag_params', '-v7.3');

fprintf('\nSaved to: %s\n', out_folder);
