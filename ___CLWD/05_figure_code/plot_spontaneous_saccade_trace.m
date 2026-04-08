function [all_fig_handles, all_stats, params] = plot_spontaneous_saccade_trace(varargin)
% plot_spontaneous_saccade_trace  v2.1
% Full-session FR + Eye Velocity trace with stationary shading, plus
% per-neuron spontaneous saccade burst modulation statistics.
%
%   Uses Chenhao's saccade onset times for burst modulation stats.
%   Eye velocity traces from saccades-edited.mat (diffHdeg/diffH).
%   Splits VOR vs Spontaneous using STtstart/STtend.
%
%   INPUTS (name-value pairs via varargin):
%     neuron_list       - cell array of neuron name strings
%     session_list      - cell array of session path strings
%     title_mode        - 'full' (stats in title) or 'simple' (neuron name only)
%     gauss_sd, fr_dt   - FR smoothing parameters
%     burst_window, baseline_window, peth_window, peth_dt, peth_smooth_win
%     min_saccades, alpha_threshold
%     doSavePNG/FIG/EPS, doShow, save_dir
%
%   OUTPUTS:
%     all_fig_handles, all_stats, params

%% =====================================================================
%  DEFAULTS
%  =====================================================================
neuron_list   = {'M086-2020-11-24-TT02_1', ...
                 'M287-2022-08-08-TT08_1', ...
                 'M284-2022-05-12-1-TT07_1'};
session_list  = {'C:\Jeff\U01\datatouse\M086\M086-2020-11-24', ...
                 'C:\Jeff\U01\datatouse\M287\M287-2022-08-08', ...
                 'C:\Jeff\U01\datatouse\M284\M284-2022-05-12-1'};

title_mode        = 'full';    % 'full' or 'simple'
gauss_sd          = 0.05;
fr_dt             = 0.001;
burst_window      = [-0.05 0.05];
baseline_window   = [-0.2 -0.1];
peth_window       = [-0.2 0.2];
peth_dt           = 0.005;
peth_smooth_win   = 5;
peth_excessBounds = 1;
min_saccades      = 3;
alpha_threshold   = 0.05;
doSavePNG         = 1;
doSaveFIG         = 1;
doSaveEPS         = 1;
doShow            = 1;
save_dir          = 'C:\Jeff\U01\datatouse\___analysis_output';
precomputed_stats = [];   % pass all_stats cell array from prior run to skip recomputation

process_varargin(varargin);

% Manual extraction for cell/struct params (process_varargin can silently fail)
for iV = 1:2:length(varargin)
    if ischar(varargin{iV})
        switch varargin{iV}
            case 'precomputed_stats'; precomputed_stats = varargin{iV+1};
            case 'neuron_list';       neuron_list = varargin{iV+1};
            case 'session_list';      session_list = varargin{iV+1};
        end
    end
end

%% =====================================================================
%  PARAMS
%  =====================================================================
params = struct();
params.function_name    = 'plot_spontaneous_saccade_trace';
params.version          = '2.1';
params.date_run         = datestr(now, 'yyyy-mm-dd HH:MM:SS');
params.neuron_list      = {neuron_list};
params.session_list     = {session_list};
params.title_mode       = title_mode;
params.gauss_sd         = gauss_sd;
params.fr_dt            = fr_dt;
params.burst_window     = burst_window;
params.baseline_window  = baseline_window;
params.peth_window      = peth_window;
params.peth_dt          = peth_dt;
params.peth_smooth_win  = peth_smooth_win;
params.min_saccades     = min_saccades;
params.alpha_threshold  = alpha_threshold;
params.saccade_source   = 'Chenhao onset times (saccade_ONSET_all_unit.mat)';

%% =====================================================================
%  SETUP
%  =====================================================================
peth_bins     = linspace(peth_window(1), peth_window(2), diff(peth_window)/peth_dt + 1);
peth_centers  = peth_bins(1:end-1) + diff(peth_bins) / 2;
nBins         = length(peth_centers);

nNeurons = length(neuron_list);
if length(session_list) ~= nNeurons
    error('neuron_list and session_list must have the same length.');
end

% Output folder
run_date_str = datestr(now, 'yyyy-mm-dd');
run_num = 1;
out_folder = fullfile(save_dir, sprintf('saccade_trace_%s_run%d', run_date_str, run_num));
while exist(out_folder, 'dir')
    run_num = run_num + 1;
    out_folder = fullfile(save_dir, sprintf('saccade_trace_%s_run%d', run_date_str, run_num));
end
if doSavePNG || doSaveFIG || doSaveEPS
    mkdir(out_folder);
    fprintf('Output folder: %s\n', out_folder);
end
params.output_folder = out_folder;

all_fig_handles = cell(nNeurons, 1);
all_stats       = cell(nNeurons, 1);

fprintf('\n============================================================\n');
fprintf('  plot_spontaneous_saccade_trace v2.1\n');
fprintf('  %d neuron(s) | saccade source: Chenhao onset times\n', nNeurons);
fprintf('============================================================\n');

%% =====================================================================
%  WORKSPACE MODE — reprint stats from prior run, no recomputation
%  =====================================================================
if ~isempty(precomputed_stats)
    fprintf('=== WORKSPACE MODE (%d entries) ===\n', length(precomputed_stats));
    all_stats = precomputed_stats;
    all_fig_handles = {};

    % Output folder for any re-saves
    run_date_str = datestr(now, 'yyyy-mm-dd');
    run_num = 1;
    out_folder = fullfile(save_dir, sprintf('saccade_trace_%s_run%d', run_date_str, run_num));
    while exist(out_folder, 'dir')
        run_num = run_num + 1;
        out_folder = fullfile(save_dir, sprintf('saccade_trace_%s_run%d', run_date_str, run_num));
    end
    params.output_folder = out_folder;

    for iN = 1:length(all_stats)
        st = all_stats{iN};
        if ~isstruct(st) || ~isfield(st, 'status'); continue; end
        if ~strcmp(st.status, 'OK'); continue; end
        fprintf('\n--- %s ---\n', st.neuron_name);
        fprintf('  Pref: %s | VOR: CCW=%d CW=%d | Spont: CCW=%d CW=%d\n', ...
            st.pref_direction, st.nCCW_move, st.nCW_move, st.nCCW_stat, st.nCW_stat);
        fprintf('  BM(VOR pref)=%.3f p=%.2e | BM(Spont pref)=%.3f p=%.2e\n', ...
            st.pref_bm_vor, st.pref_p_vor, st.pref_bm_spont, st.pref_p_spont);
    end
    fprintf('\nWorkspace mode: stats printed. Pass no precomputed_stats to recompute + plot.\n');
    return
end

%% =====================================================================
%  MAIN LOOP
%  =====================================================================
for iN = 1:nNeurons
    neuron_name  = neuron_list{iN};
    session_path = session_list{iN};
    fprintf('\n--- %d/%d: %s ---\n', iN, nNeurons, neuron_name);

    stat_struct = struct();
    stat_struct.neuron_name  = neuron_name;
    stat_struct.session_path = session_path;

    try
        pushdir(session_path);
        SSN = HD_GetSSN;
        stat_struct.SSN = SSN;

        %% Load spikes + compute FR via MakeQfromS
        [S_all, spike_cell_names] = LoadSpikesQuick();
        cell_idx = find_cell_index(spike_cell_names, neuron_name, SSN);
        if isempty(cell_idx)
            fprintf('  ERROR: cell not found\n');
            stat_struct.status = 'ERROR: cell not found';
            all_stats{iN} = stat_struct; popdir; continue
        end
        myCell = SelectTS([], S_all, cell_idx);

        cfg_Q = []; cfg_Q.dt = fr_dt; cfg_Q.gausswin_sd = gauss_sd; cfg_Q.smooth = 'gauss';
        Q = MakeQfromS(cfg_Q, myCell);
        FR_trace = Q.data ./ cfg_Q.dt;
        FR_tvec  = Q.tvec;

        %% Load eye velocity (from saccades-edited.mat — still needed for traces)
        sacc_file = strcat(SSN, '-saccades-edited.mat');
        if ~exist(sacc_file, 'file')
            stat_struct.status = 'ERROR: no saccades file';
            all_stats{iN} = stat_struct; popdir; continue
        end
        eye_fields = load(sacc_file, 'diffHdeg', 'diffH');
        if isfield(eye_fields, 'diffHdeg')
            EV = eye_fields.diffHdeg;
        elseif isfield(eye_fields, 'diffH')
            EV = eye_fields.diffH;
            fprintf('  WARNING: using diffH (pixels/s)\n');
        else
            stat_struct.status = 'ERROR: no EV data';
            all_stats{iN} = stat_struct; popdir; continue
        end

        %% Load Chenhao's saccade onset times
        [temporal_all, nasal_all] = load_chenhao_saccades(SSN);
        if isempty(temporal_all) && isempty(nasal_all)
            fprintf('  WARNING: no Chenhao saccade times for %s\n', SSN);
            stat_struct.status = 'ERROR: no Chenhao saccades';
            all_stats{iN} = stat_struct; popdir; continue
        end

        %% Load stationary times
        stat_file = strcat(SSN, '-AHV_StationaryTimes.mat');
        if ~exist(stat_file, 'file')
            stat_struct.status = 'ERROR: no StationaryTimes';
            all_stats{iN} = stat_struct; popdir; continue
        end
        stat_data = load(stat_file, 'STtstart', 'STtend');
        STtstart = stat_data.STtstart(:);
        STtend   = stat_data.STtend(:);
        total_stat_dur = sum(STtend - STtstart);

        %% Split VOR vs stationary (matching v1.3)
        t_stat_mask = false(size(temporal_all));
        n_stat_mask = false(size(nasal_all));
        for iST = 1:length(STtstart)
            t_stat_mask = t_stat_mask | (temporal_all > STtstart(iST) & temporal_all < STtend(iST));
            n_stat_mask = n_stat_mask | (nasal_all > STtstart(iST) & nasal_all < STtend(iST));
        end
        temporal_stat = temporal_all(t_stat_mask);
        nasal_stat    = nasal_all(n_stat_mask);
        temporal_move = temporal_all(~t_stat_mask);
        nasal_move    = nasal_all(~n_stat_mask);

        stat_struct.nCCW_stat = length(temporal_stat);
        stat_struct.nCW_stat  = length(nasal_stat);
        stat_struct.nCCW_move = length(temporal_move);
        stat_struct.nCW_move  = length(nasal_move);
        stat_struct.stat_dur  = total_stat_dur;

        %% Compute burst modulation (matching v1.3 exactly)
        cond_labels = {'move_CCW', 'move_CW', 'stat_CCW', 'stat_CW'};
        cond_sacc   = {temporal_move, nasal_move, temporal_stat, nasal_stat};
        cond_n      = [length(temporal_move), length(nasal_move), length(temporal_stat), length(nasal_stat)];

        for iCond = 1:4
            if cond_n(iCond) < 1
                stat_struct.(cond_labels{iCond}) = struct('fr_smooth', nan(1,nBins), ...
                    'burst_mod', NaN, 'p_value', NaN, 'is_sig', false, ...
                    'mean_burst_fr', NaN, 'mean_baseline_fr', NaN, 'nEvents', 0);
                continue
            end

            % Raster + PETH
            spk_rel_all = [];
            for iE = 1:length(cond_sacc{iCond})
                S0 = restrict(myCell, cond_sacc{iCond}(iE)+peth_window(1)-peth_excessBounds, ...
                                      cond_sacc{iCond}(iE)+peth_window(2)+peth_excessBounds);
                S0 = restrict(S0, cond_sacc{iCond}(iE)+peth_window(1), ...
                                  cond_sacc{iCond}(iE)+peth_window(2));
                if ~isempty(S0.t{1})
                    spk_rel_all = [spk_rel_all; S0.t{1} - cond_sacc{iCond}(iE)]; %#ok<AGROW>
                end
            end
            m_hist = histcounts(spk_rel_all, peth_bins);
            fr_peth = m_hist / peth_dt / max(cond_n(iCond), 1);
            fr_smooth = smoothdata(movmean(fr_peth, peth_smooth_win), 'gaussian', peth_smooth_win);

            % Per-trial burst vs baseline
            trial_burst_fr = nan(cond_n(iCond), 1);
            trial_baseline_fr = nan(cond_n(iCond), 1);
            burst_dur = diff(burst_window);
            base_dur  = diff(baseline_window);
            for iE = 1:cond_n(iCond)
                S_b  = restrict(myCell, cond_sacc{iCond}(iE)+burst_window(1), cond_sacc{iCond}(iE)+burst_window(2));
                trial_burst_fr(iE) = length(S_b.t{1}) / burst_dur;
                S_bl = restrict(myCell, cond_sacc{iCond}(iE)+baseline_window(1), cond_sacc{iCond}(iE)+baseline_window(2));
                trial_baseline_fr(iE) = length(S_bl.t{1}) / base_dur;
            end
            mean_burst = mean(trial_burst_fr, 'omitnan');
            mean_base  = mean(trial_baseline_fr, 'omitnan');
            if mean_base > 0; bm_val = (mean_burst - mean_base) / mean_base; else; bm_val = NaN; end
            valid_tr = ~isnan(trial_burst_fr) & ~isnan(trial_baseline_fr);
            if sum(valid_tr) >= 3; [p_val, ~] = signrank(trial_burst_fr(valid_tr), trial_baseline_fr(valid_tr));
            else; p_val = NaN; end

            stat_struct.(cond_labels{iCond}) = struct('fr_smooth', fr_smooth, 'burst_mod', bm_val, ...
                'p_value', p_val, 'is_sig', p_val < alpha_threshold, ...
                'mean_burst_fr', mean_burst, 'mean_baseline_fr', mean_base, 'nEvents', cond_n(iCond));
        end

        %% Preferred direction
        mod_CCW_vor = stat_struct.move_CCW.burst_mod;
        mod_CW_vor  = stat_struct.move_CW.burst_mod;
        if isnan(mod_CCW_vor) && isnan(mod_CW_vor); pref_dir = 'none';
        elseif isnan(mod_CW_vor) || mod_CCW_vor >= mod_CW_vor; pref_dir = 'CCW';
        else; pref_dir = 'CW'; end
        stat_struct.pref_direction = pref_dir;

        if strcmp(pref_dir, 'CCW'); ps = 'CCW'; nps = 'CW';
        else; ps = 'CW'; nps = 'CCW'; end
        stat_struct.pref_bm_vor   = stat_struct.(['move_' ps]).burst_mod;
        stat_struct.pref_bm_spont = stat_struct.(['stat_' ps]).burst_mod;
        stat_struct.pref_p_vor    = stat_struct.(['move_' ps]).p_value;
        stat_struct.pref_p_spont  = stat_struct.(['stat_' ps]).p_value;

        %% Print stats
        fprintf('\n  Session:    %s\n', SSN);
        fprintf('  Pref dir:   %s\n', pref_dir);
        fprintf('  Stationary: %.1f s total across %d periods\n', total_stat_dur, length(STtstart));
        fprintf('  Saccade Counts:\n');
        fprintf('    VOR:          CCW=%d  CW=%d\n', cond_n(1), cond_n(2));
        fprintf('    Spontaneous:  CCW=%d  CW=%d\n', cond_n(3), cond_n(4));
        fprintf('  Burst Modulation (pref=%s):\n', pref_dir);
        fprintf('    VOR pref:     BM=%.3f  p=%.4e  %s\n', stat_struct.pref_bm_vor, ...
            stat_struct.pref_p_vor, sig_str(stat_struct.pref_p_vor, alpha_threshold));
        fprintf('    Spont pref:   BM=%.3f  p=%.4e  %s\n', stat_struct.pref_bm_spont, ...
            stat_struct.pref_p_spont, sig_str(stat_struct.pref_p_spont, alpha_threshold));
        fprintf('  All Conditions:\n');
        fprintf('    %-12s  %6s  %8s  %10s  %5s\n', 'Condition', 'N', 'BurstMod', 'p-value', 'Sig');
        fprintf('    %s\n', repmat('-', 1, 50));
        for iCond = 1:4
            cd_s = stat_struct.(cond_labels{iCond});
            fprintf('    %-12s  %6d  %8.3f  %10.4e  %5s\n', cond_labels{iCond}, ...
                cd_s.nEvents, cd_s.burst_mod, cd_s.p_value, sig_str(cd_s.p_value, alpha_threshold));
        end

        stat_struct.status = 'OK';

        %% =============================================================
        %  PLOT
        %  =============================================================
        if doShow; fh = figure('Position', [50 50 1900 550]);
        else; fh = figure('Position', [50 50 1900 550], 'Visible', 'off'); end

        fh.UserData.date_created        = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        fh.UserData.generating_function = 'plot_spontaneous_saccade_trace v2.1';
        fh.UserData.neuron_name         = neuron_name;
        hold on

        % FR on left y-axis
        yyaxis left
        h_FR = plot(FR_tvec, FR_trace, '-', 'Color', [0.20 0.50 0.80], 'LineWidth', 1.2);
        ylabel('Firing Rate (Hz)', 'FontSize', 22, 'FontWeight', 'bold')
        set(gca, 'FontSize', 18)
        FR_color = h_FR.Color;

        % EV on right y-axis
        yyaxis right
        h_EV = plot(EV.tvec, EV.data, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 0.6);
        ylabel('Eye Velocity (deg/s)', 'FontSize', 22, 'FontWeight', 'bold')
        EV_color = h_EV.Color;

        xlabel('Time (s)', 'FontSize', 22, 'FontWeight', 'bold')

        % Title based on title_mode
        if strcmp(title_mode, 'simple')
            title(neuron_name, 'FontSize', 24, 'Interpreter', 'none', 'FontWeight', 'bold')
        else
            title_str = sprintf('%s  |  Pref: %s  |  BM(VOR)=%.2f  BM(Spont)=%.2f', ...
                neuron_name, pref_dir, stat_struct.pref_bm_vor, stat_struct.pref_bm_spont);
            title(title_str, 'FontSize', 22, 'Interpreter', 'none', 'FontWeight', 'bold')
        end

        % Stationary shading
        yyaxis left
        ax_limits = axis;
        y_lo = ax_limits(3); y_hi = ax_limits(4);
        for iST = 1:length(STtstart)
            t0 = STtstart(iST); t1 = STtend(iST);
            if ~any(isnan([t0 t1 y_lo y_hi]))
                fill([t0 t1 t1 t0], [y_lo y_lo y_hi y_hi], [0.85 0.85 0.85], ...
                    'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
        end

        % Replot on top
        yyaxis left
        plot(FR_tvec, FR_trace, '-', 'Color', FR_color, 'LineWidth', 1.2);
        yyaxis right
        plot(EV.tvec, EV.data, '-', 'Color', EV_color, 'LineWidth', 0.6);

        gca_handle = gca; gca_handle.Box = 'off';

        % Legend
        yyaxis left
        h_leg = fill(nan, nan, [0.85 0.85 0.85], 'EdgeColor', 'none');
        legend(h_leg, 'Stationary', 'FontSize', 16, 'Location', 'best', 'Box', 'off');

        %% Save
        if doSavePNG || doSaveFIG || doSaveEPS
            safe_name = strrep(neuron_name, '-', '_');
            base_fname = sprintf('session_FR_EV_%s', safe_name);
            save_fig_trio(fh, fullfile(out_folder, base_fname), doSavePNG, doSaveFIG, doSaveEPS, doShow);
        end

        all_fig_handles{iN} = fh;
        popdir;

    catch ME
        stat_struct.status = sprintf('ERROR: %s', ME.message);
        fprintf('  ERROR: %s\n', ME.message);
        try popdir; catch; end
    end
    all_stats{iN} = stat_struct;
end

%% Save stats
if doSavePNG || doSaveFIG || doSaveEPS
    save(fullfile(out_folder, 'trace_stats.mat'), 'all_stats', 'params');
    fprintf('\nSaved to: %s\n', out_folder);
end
fprintf('\nDone. %d neuron(s) processed.\n', nNeurons);
end

%% =====================================================================
%  LOCAL HELPERS
%  =====================================================================
function save_fig_trio(fh, fp, doP, doF, doE, doS)
    if doP; print(fh, [fp '.png'], '-dpng', '-r300'); end
    if doF; savefig(fh, [fp '.fig']); end
    if doE; print(fh, [fp '.eps'], '-depsc2'); end
    if ~doS; close(fh); end
end

function cell_idx = find_cell_index(names, neuron, SSN)
    cell_idx = [];
    for cc = 1:length(names)
        if strcmp(names{cc}, neuron); cell_idx = cc; return; end
    end
    for cc = 1:length(names)
        if contains(names{cc}, strrep(neuron, [SSN '-'], '')); cell_idx = cc; return; end
    end
end

function s = sig_str(p, alpha)
    if isnan(p); s = ''; return; end
    if p < alpha/10; s = '***'; elseif p < alpha; s = '*'; else; s = 'n.s.'; end
end
