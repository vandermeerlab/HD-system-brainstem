%% run_optotag_eye_shuffle.m  — v1.1
%
% Runs compute_optotag_eye_shuffle on all 7 optotagged sessions using
% BOTH stationary-restricted shuffle methods:
%   1) 'random_stationary'   — random placement within stationary periods
%   2) 'circular_stationary' — circular shift within concatenated stationary timeline
%
% Results for both methods are saved to workspace and disk.
%
% TWO MODES:
%   1) FULL RUN — set replot_only = 0 (default)
%   2) REPLOT ONLY — set replot_only = 1, requires results in workspace
%
% WORKSPACE VARIABLES AFTER FULL RUN:
%   eye_shuff_random   — cell array (7x1) of shuffle_out structs (method 1)
%   eye_shuff_circular — cell array (7x1) of shuffle_out structs (method 2)
%   eye_shuff_params   — cell array (7x1) of params structs (from method 1)
%   session_status     — cell array (7x1) of status strings
%   opto_sessions_list, opto_neurons_list, runner_params
%
% REPLOT:
%   replot_only = 1; run_optotag_eye_shuffle
%
% LOAD FROM DISK:
%   load('...\eye_shuffle_results.mat');
%   replot_only = 1; run_optotag_eye_shuffle
%
% Version History:
%   2026-04-07  v1.1  Both shuffle methods. Stationary-period restriction.
%                     500 shuffles. Workspace replot mode.
%   2026-04-07  v1.0  Initial version. Unrestricted shuffle only.

%% =====================================================================
%  CONFIG
%  =====================================================================
if ~exist('replot_only', 'var'); replot_only = 0; end

nShuffles_cfg = 500;

data_root = 'C:\Jeff\U01\datatouse';
save_dir  = 'C:\Jeff\U01\datatouse\___analysis_output';

opto_sessions_list = {'M271-2021-08-30', 'M282-2022-02-10', 'M282-2022-02-15', ...
    'M391-2022-12-07-1', 'M391-2022-12-07-2', 'M402-2022-11-02-2', ...
    'M444-2023-08-18'};
opto_neurons_list = {'M271-2021-08-30-TT02_1', 'M282-2022-02-10-TT03_2', ...
    'M282-2022-02-15-TT07_1', 'M391-2022-12-07-1-TT01_1', ...
    'M391-2022-12-07-2-TT02_1', 'M402-2022-11-02-2-TT02_3', ...
    'M444-2023-08-18-TT07_1'};

nSessions = length(opto_sessions_list);

%% =====================================================================
%  OUTPUT FOLDER
%  =====================================================================
run_date_str = datestr(now, 'yyyy-mm-dd');
run_time_str = datestr(now, 'hh-MM-SS_AM');
run_folder_name = sprintf('optotag_eye_shuffle_%s_%s', run_date_str, run_time_str);
out_folder = fullfile(save_dir, run_folder_name);

%% =====================================================================
%  FULL RUN
%  =====================================================================
if ~replot_only
    mkdir(out_folder);
    fprintf('Output folder: %s\n', out_folder);

    eye_shuff_random   = cell(nSessions, 1);
    eye_shuff_circular = cell(nSessions, 1);
    eye_shuff_params   = cell(nSessions, 1);
    session_status     = cell(nSessions, 1);

    for iS = 1:nSessions
        SSN = opto_sessions_list{iS};
        mouse_id = SSN(1:4);
        sess_path = fullfile(data_root, mouse_id, SSN);
        fprintf('\n=== %d/%d: %s ===\n', iS, nSessions, SSN);

        % Laser duration per session
        if contains(SSN, 'M271')
            this_laser_dur = 0.5;
        else
            this_laser_dur = 1.0;
        end

        try
            pushdir(sess_path);

            % --- Load laser times ---
            [~, ~, laser_on_times, ~, n_pulses] = SortBrainstemEventLabels3;
            if isnan(laser_on_times)
                session_status{iS} = 'SKIP: no laser events';
                fprintf('  No laser events — skipping\n');
                popdir; continue
            end

            % --- Get session duration ---
            events_ts = LoadEvents([]);
            start_label_idx = [];
            for iL = 1:length(events_ts.label)
                if strcmp(events_ts.label{iL}, 'Starting Recording')
                    start_label_idx = iL; break
                end
            end
            stop_label_idx = [];
            for iL = 1:length(events_ts.label)
                if strcmp(events_ts.label{iL}, 'Stopping Recording')
                    stop_label_idx = iL; break
                end
            end
            start_ts = events_ts.t{start_label_idx}(1);
            if ~isempty(stop_label_idx)
                session_dur = events_ts.t{stop_label_idx}(end) - start_ts;
            else
                session_dur = max(laser_on_times) + 60;
                fprintf('  WARNING: no stop time, using fallback session duration\n');
            end

            % --- Load stationary times ---
            stat_file = strcat(SSN, '-AHV_StationaryTimes.mat');
            if ~exist(stat_file, 'file')
                session_status{iS} = 'SKIP: no StationaryTimes.mat';
                fprintf('  No StationaryTimes.mat — skipping\n');
                popdir; continue
            end
            st = load(stat_file, 'STtstart', 'STtend');
            stat_intervals = [st.STtstart(:), st.STtend(:)];
            fprintf('  %d stationary intervals (%.0f s total)\n', ...
                size(stat_intervals, 1), sum(stat_intervals(:,2) - stat_intervals(:,1)));

            % --- Load eye position data ---
            sacc_file = strcat(SSN, '-saccades-edited.mat');
            if ~exist(sacc_file, 'file')
                session_status{iS} = 'SKIP: no saccades-edited.mat';
                fprintf('  No saccades-edited.mat — skipping\n');
                popdir; continue
            end

            eye_vars = load(sacc_file, 'tsdH', 'k');
            if ~isfield(eye_vars, 'tsdH')
                session_status{iS} = 'SKIP: no tsdH';
                fprintf('  No tsdH — skipping\n');
                popdir; continue
            end

            ep_tvec_loaded = eye_vars.tsdH.tvec(:);
            ep_data_loaded = eye_vars.tsdH.data(:);

            % Convert pixels to degrees
            if isfield(eye_vars, 'k')
                cal_k = eye_vars.k;
            else
                cal_k = 1;
                fprintf('  WARNING: no calibration factor k, eye pos in pixels\n');
            end
            ep_data_degrees = ep_data_loaded * cal_k;

            % tsdH.tvec is ALREADY session-relative (processPupilData2
            % subtracts starttime). Do NOT subtract again.

            fprintf('  %d laser pulses, %.0f s session, %d eye samples\n', ...
                n_pulses, session_dur, length(ep_tvec_loaded));

            % --- Method 1: random_stationary ---
            fprintf('\n  --- Method 1: random_stationary ---\n');
            cfg_rand = [];
            cfg_rand.nShuffles      = nShuffles_cfg;
            cfg_rand.laser_dur      = this_laser_dur;
            cfg_rand.shuffle_method = 'random_stationary';
            cfg_rand.stat_intervals = stat_intervals;
            cfg_rand.verbose        = 1;

            [eye_shuff_random{iS}, eye_shuff_params{iS}] = compute_optotag_eye_shuffle( ...
                ep_tvec_loaded, ep_data_degrees, laser_on_times, session_dur, cfg_rand);

            % --- Method 2: circular_stationary ---
            fprintf('\n  --- Method 2: circular_stationary ---\n');
            cfg_circ = [];
            cfg_circ.nShuffles      = nShuffles_cfg;
            cfg_circ.laser_dur      = this_laser_dur;
            cfg_circ.shuffle_method = 'circular_stationary';
            cfg_circ.stat_intervals = stat_intervals;
            cfg_circ.verbose        = 1;

            [eye_shuff_circular{iS}, ~] = compute_optotag_eye_shuffle( ...
                ep_tvec_loaded, ep_data_degrees, laser_on_times, session_dur, cfg_circ);

            session_status{iS} = 'OK';
            fprintf('\n  RANDOM:   pctile=%.4f  peak=%.3f deg\n', ...
                eye_shuff_random{iS}.percentile_overall, eye_shuff_random{iS}.peak_disp_real);
            fprintf('  CIRCULAR: pctile=%.4f  peak=%.3f deg\n', ...
                eye_shuff_circular{iS}.percentile_overall, eye_shuff_circular{iS}.peak_disp_real);

            popdir;
        catch ME
            session_status{iS} = sprintf('ERROR: %s', ME.message);
            fprintf('  ERROR: %s\n', ME.message);
            try popdir; catch; end
        end
    end

    %% Save to disk
    runner_params = struct();
    runner_params.function_name = 'run_optotag_eye_shuffle';
    runner_params.version       = '1.1';
    runner_params.date_run      = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    runner_params.nShuffles     = nShuffles_cfg;
    runner_params.methods_run   = {'random_stationary', 'circular_stationary'};
    runner_params.opto_sessions = {opto_sessions_list};
    runner_params.opto_neurons  = {opto_neurons_list};

    save(fullfile(out_folder, 'eye_shuffle_results.mat'), ...
        'eye_shuff_random', 'eye_shuff_circular', 'eye_shuff_params', ...
        'session_status', 'opto_sessions_list', 'opto_neurons_list', ...
        'runner_params', '-v7.3');
    fprintf('\nData saved to: %s\n', out_folder);

else
    %% REPLOT MODE
    if ~exist('eye_shuff_random', 'var') || ~exist('eye_shuff_circular', 'var')
        error('replot_only=1 but eye_shuff_random/circular not in workspace.');
    end
    if ~exist('session_status', 'var')
        error('replot_only=1 but session_status not in workspace.');
    end
    fprintf('=== REPLOT MODE — using workspace data ===\n');
    mkdir(out_folder);
    fprintf('Replot output folder: %s\n', out_folder);
end

%% =====================================================================
%  SUMMARY TABLE (both modes)
%  =====================================================================
fprintf('\n============================================================\n');
fprintf('LASER EYE DISPLACEMENT SHUFFLE — BOTH METHODS\n');
fprintf('============================================================\n');
fprintf('%-25s  %12s  %12s  %12s  %12s\n', ...
    'Session', 'Rand pctile', 'Rand peak', 'Circ pctile', 'Circ peak');
fprintf('%-25s  %12s  %12s  %12s  %12s\n', ...
    repmat('-',1,25), repmat('-',1,12), repmat('-',1,12), repmat('-',1,12), repmat('-',1,12));

for iS = 1:nSessions
    if strcmp(session_status{iS}, 'OK')
        sr = eye_shuff_random{iS};
        sc = eye_shuff_circular{iS};
        sig_r = ''; sig_c = '';
        if sr.percentile_overall >= 0.95; sig_r = ' *'; end
        if sr.percentile_overall >= 0.99; sig_r = ' **'; end
        if sr.percentile_overall >= 0.999; sig_r = '***'; end
        if sc.percentile_overall >= 0.95; sig_c = ' *'; end
        if sc.percentile_overall >= 0.99; sig_c = ' **'; end
        if sc.percentile_overall >= 0.999; sig_c = '***'; end
        fprintf('%-25s  %9.4f%3s  %9.3f     %9.4f%3s  %9.3f\n', ...
            opto_sessions_list{iS}, ...
            sr.percentile_overall, sig_r, sr.peak_disp_real, ...
            sc.percentile_overall, sig_c, sc.peak_disp_real);
    else
        fprintf('%-25s  %s\n', opto_sessions_list{iS}, session_status{iS});
    end
end
fprintf('\n* p<.05   ** p<.01   *** p<.001\n');

%% =====================================================================
%  FIGURE 1: RANDOM STATIONARY METHOD
%  =====================================================================
plot_summary_figure(eye_shuff_random, eye_shuff_params, session_status, ...
    opto_sessions_list, nSessions, out_folder, ...
    'Random Stationary Shuffle', 'eye_shuffle_random_stationary');

%% =====================================================================
%  FIGURE 2: CIRCULAR STATIONARY METHOD
%  =====================================================================
plot_summary_figure(eye_shuff_circular, eye_shuff_params, session_status, ...
    opto_sessions_list, nSessions, out_folder, ...
    'Circular Stationary Shuffle', 'eye_shuffle_circular_stationary');

fprintf('\nDone. Both figures saved to: %s\n', out_folder);

%% =====================================================================
%  LOCAL FUNCTION: summary figure
%  =====================================================================
function plot_summary_figure(shuff_results, shuff_params_cells, sess_status, ...
        sess_names, nSess, folder, fig_title_str, file_prefix)

    fig_h = figure('Position', [50 50 1800 900], 'Visible', 'on');
    fig_h.UserData.date_created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fig_h.UserData.generating_function = 'run_optotag_eye_shuffle v1.1';

    n_ok = sum(strcmp(sess_status, 'OK'));
    if n_ok == 0
        text(0.5, 0.5, 'No valid sessions', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 24);
    else
        plot_counter = 0;
        for iS = 1:nSess
            if ~strcmp(sess_status{iS}, 'OK'); continue; end
            plot_counter = plot_counter + 1;
            sr = shuff_results{iS};

            % Laser duration
            if ~isempty(shuff_params_cells{iS})
                ld = shuff_params_cells{iS}.laser_dur;
            else
                ld = 1.0;
            end

            % --- Top row: displacement traces ---
            ax_t = subplot(2, n_ok, plot_counter); hold(ax_t, 'on');

            % Shuffle mean +/- SEM
            fx = [sr.eye_grid fliplr(sr.eye_grid)];
            fy = [sr.shuff_mean + sr.shuff_sem, fliplr(sr.shuff_mean - sr.shuff_sem)];
            if ~any(isnan(fy))
                fill(ax_t, fx, fy, [0.7 0.7 0.7], ...
                    'EdgeColor', 'none', 'FaceAlpha', 0.4);
            end
            plot(ax_t, sr.eye_grid, sr.shuff_mean, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
            plot(ax_t, sr.eye_grid, sr.disp_real_smooth, '-', 'Color', [0.60 0.15 0.60], 'LineWidth', 3);

            % Laser shading
            yl_ax = ylim(ax_t);
            fill(ax_t, [0 ld ld 0], [yl_ax(1) yl_ax(1) yl_ax(2) yl_ax(2)], ...
                [0.4 0.7 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.15);

            % Sig bars
            sig_any = sr.sig_001 | sr.sig_05;
            if any(sig_any)
                yl_bar = yl_ax(2) * 0.95;
                sig_x = sr.eye_grid(sig_any);
                plot(ax_t, sig_x, repmat(yl_bar, 1, length(sig_x)), 's', ...
                    'MarkerSize', 4, 'MarkerFaceColor', [0 0.6 0], 'MarkerEdgeColor', 'none');
            end

            xline(ax_t, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1);
            xlim(ax_t, [-1 2]); box(ax_t, 'off');
            title(ax_t, sprintf('%s\np=%.3f', sess_names{iS}, ...
                1 - sr.percentile_overall), 'FontSize', 14, 'Interpreter', 'none');
            if plot_counter == 1; ylabel(ax_t, '|Eye Disp| (deg)', 'FontSize', 18); end
            set(ax_t, 'FontSize', 14);

            % --- Bottom row: histogram ---
            ax_h = subplot(2, n_ok, n_ok + plot_counter); hold(ax_h, 'on');
            histogram(ax_h, sr.peak_disp_shuff, 30, 'FaceColor', [0.7 0.7 0.7], ...
                'EdgeColor', 'none', 'FaceAlpha', 0.7);
            xline(ax_h, sr.peak_disp_real, '-', 'Color', [0.60 0.15 0.60], 'LineWidth', 3);
            xlabel(ax_h, 'Peak |Disp| (deg)', 'FontSize', 16);
            if plot_counter == 1; ylabel(ax_h, 'Count', 'FontSize', 18); end
            title(ax_h, sprintf('%.1f%% ile', sr.percentile_overall * 100), 'FontSize', 14);
            set(ax_h, 'FontSize', 14); box(ax_h, 'off');
        end
    end
    sgtitle(sprintf('Laser Eye Displacement — %s', fig_title_str), ...
        'FontSize', 24, 'FontWeight', 'bold');

    print(fig_h, fullfile(folder, [file_prefix '.png']), '-dpng', '-r300');
    savefig(fig_h, fullfile(folder, [file_prefix '.fig']));
    fprintf('Saved: %s.png\n', file_prefix);
end
