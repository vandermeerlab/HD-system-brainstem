function [batch_summary, all_results, params] = analyze_stationary_saccade_bursts_ver_1_2(varargin)
% analyze_stationary_saccade_bursts_ver_1_2  Population analysis of
% saccade-evoked bursting: VOR (moving) vs. Spontaneous (stationary).
%
%   Compares saccade-aligned neural activity during VOR-evoked saccades
%   (platform moving) vs. spontaneous saccades (platform stationary).
%   For optotagged sessions, laser bout periods are excluded from
%   stationary saccade counts to prevent contamination.
%
%   Generates two population figures:
%     1. Z-scored population PETH (preferred vs non-preferred x VOR vs spontaneous)
%     2. Raw FR population PETH (same layout, no normalization)
%
%   Also outputs: best cells list, saccade imbalance flags, significance
%   tests (preferred vs non-preferred, VOR vs spontaneous).
%
%   INPUTS (name-value pairs via varargin):
%       cell_list           - cell array of full t-file paths
%       precomputed_results - cell array from prior run (workspace mode)
%       pref_direction_list - containers.Map: neuron -> 'CCW'|'CW'
%       min_saccades        - both dirs must meet this (default: 8)
%       burst_window        - [start end] s (default: [-0.05 0.05])
%       baseline_window     - [start end] s (default: [-0.2 -0.1])
%       laser_bout_buffer   - s, exclusion buffer around laser bouts (default: 1.0)
%       imbalance_ratio     - flag sessions with ratio > this (default: 3.0)
%       n_best_cells        - how many top cells to report (default: 15)
%       doMakeFigs, doPopulationFig, doSavePNG, doSaveFIG, doShow
%
%   OUTPUTS:
%       batch_summary   - MATLAB table
%       all_results     - cell array of per-neuron structs
%       params          - struct of parameter values
%
%   USAGE:
%       % Full run:
%       [s, r, p] = analyze_stationary_saccade_bursts_ver_1_2();
%
%       % Workspace mode from prior run:
%       [~, ~, p] = analyze_stationary_saccade_bursts_ver_1_2(...
%           'precomputed_results', all_results, 'min_saccades', 10, 'doShow', 1);
%
%   DEPENDENCIES:
%       LoadSpikesQuick, SelectTS, restrict, pushdir, popdir, HD_GetSSN,
%       process_varargin, SortBrainstemEventLabels3

%% ========================================================================
%  DEFAULTS
%  ========================================================================
cell_list           = {};
tfilelist_path      = 'C:\Jeff\U01\datatouse\NPH all tfilelist.mat';
save_root           = 'C:\Jeff\U01\datatouse\___analysis_output';
precomputed_results = {};
pref_direction_list = [];

% Analysis windows
burst_window        = [-0.05 0.05];
baseline_window     = [-0.2 -0.1];
alpha_threshold     = 0.05;

% PETH parameters
peth_window         = [-0.2 0.2];
peth_dt             = 0.005;
peth_smooth_win     = 5;
peth_excessBounds   = 1;

% Thresholds
min_saccades        = 8;
imbalance_ratio     = 3.0;
n_best_cells        = 15;

% Laser bout exclusion
laser_bout_buffer   = 1.0;     % s buffer before first pulse and after last pulse offset

% Known optotagged sessions (for laser exclusion)
opto_sessions = {'M271-2021-08-30', 'M282-2022-02-10', 'M282-2022-02-15', ...
    'M391-2022-12-07-1', 'M391-2022-12-07-2', 'M402-2022-11-02-2', ...
    'M444-2023-08-18'};

% Figure parameters
FontSize            = 18;
TickFontSize        = 14;
TitleFontSize       = 16;
MarkerSize_raster   = 3;
LineWidth_FR        = 4;
color_vor_pref      = [0.20 0.50 0.80];
color_vor_nonpref   = [0.90 0.40 0.15];
color_spont_pref    = [0.08 0.25 0.50];
color_spont_nonpref = [0.60 0.22 0.05];
rasterColor         = [0.25 0.25 0.25];
figWidth            = 1200;
figHeight           = 900;

% Output control
doMakeFigs      = 1;
doPopulationFig = 1;
doSavePNG       = 1;
doSaveFIG       = 1;
doShow          = 0;
verbose         = 1;

%% Process varargin
process_varargin(varargin);

%% ========================================================================
%  PARAMS
%  ========================================================================
params = struct();
params.function_name        = 'analyze_stationary_saccade_bursts';
params.version              = '1.2';
params.date_run             = datestr(now, 'yyyy-mm-dd HH:MM:SS');
params.burst_window         = burst_window;
params.baseline_window      = baseline_window;
params.alpha_threshold      = alpha_threshold;
params.peth_window          = peth_window;
params.peth_dt              = peth_dt;
params.peth_smooth_win      = peth_smooth_win;
params.min_saccades         = min_saccades;
params.laser_bout_buffer    = laser_bout_buffer;
params.imbalance_ratio      = imbalance_ratio;
params.comparison           = 'VOR (moving) vs Spontaneous (stationary)';

peth_bins    = linspace(peth_window(1), peth_window(2), diff(peth_window)/peth_dt + 1);
peth_centers = peth_bins(1:end-1) + diff(peth_bins) / 2;
nBins        = length(peth_centers);
baseline_bins_mask = peth_centers >= baseline_window(1) & peth_centers <= baseline_window(2);
burst_bins_mask    = peth_centers >= burst_window(1) & peth_centers <= burst_window(2);

%% ========================================================================
%  WORKSPACE MODE
%  ========================================================================
if ~isempty(precomputed_results)
    if verbose; fprintf('=== WORKSPACE MODE (%d entries) ===\n', length(precomputed_results)); end

    % Create output folder
    date_str = datestr(now, 'yyyy-mm-dd');
    run_num = 1;
    out_folder = fullfile(save_root, sprintf('stationary_burst_analysis_%s_run%d', date_str, run_num));
    while exist(out_folder, 'dir'); run_num = run_num + 1;
        out_folder = fullfile(save_root, sprintf('stationary_burst_analysis_%s_run%d', date_str, run_num)); end
    mkdir(out_folder);

    % Detect format
    ws_format = detect_results_format(precomputed_results);
    if verbose; fprintf('  Format: %s\n', ws_format); end

    % Get peth_centers
    for ii = 1:length(precomputed_results)
        if isfield(precomputed_results{ii}, 'peth_centers') && ~isempty(precomputed_results{ii}.peth_centers)
            peth_centers = precomputed_results{ii}.peth_centers; nBins = length(peth_centers);
            baseline_bins_mask = peth_centers >= baseline_window(1) & peth_centers <= baseline_window(2);
            burst_bins_mask = peth_centers >= burst_window(1) & peth_centers <= burst_window(2);
            break
        end
    end

    % Collect qualifying neurons
    ws_pref_vor = []; ws_pref_spont = []; ws_nonpref_vor = []; ws_nonpref_spont = [];
    ws_pref_vor_raw = []; ws_pref_spont_raw = []; ws_nonpref_vor_raw = []; ws_nonpref_spont_raw = [];
    ws_names = {}; ws_sessions = {}; ws_count = 0;

    for iN = 1:length(precomputed_results)
        nd_ws = precomputed_results{iN};
        if ~isstruct(nd_ws); continue; end

        [nT_s, nN_s, fr_T_vor, fr_N_vor, fr_T_spont, fr_N_spont, valid_ws] = ...
            extract_fr_traces(nd_ws, ws_format);
        if ~valid_ws; continue; end
        if nT_s < min_saccades || nN_s < min_saccades; continue; end
        if all(isnan(fr_T_vor)) || all(isnan(fr_N_vor)); continue; end
        if all(isnan(fr_T_spont)) || all(isnan(fr_N_spont)); continue; end

        % Determine preferred direction
        mod_T = compute_burst_mod_from_trace(fr_T_vor, baseline_bins_mask, burst_bins_mask);
        mod_N = compute_burst_mod_from_trace(fr_N_vor, baseline_bins_mask, burst_bins_mask);
        if isnan(mod_T) && isnan(mod_N); continue; end

        if ~isempty(pref_direction_list) && isa(pref_direction_list, 'containers.Map') ...
                && pref_direction_list.isKey(nd_ws.neuron_name)
            if strcmp(pref_direction_list(nd_ws.neuron_name), 'CW')
                mod_T = 0; mod_N = 1;
            else; mod_T = 1; mod_N = 0;
            end
        end

        ws_count = ws_count + 1;
        ws_names{ws_count} = nd_ws.neuron_name;
        if isfield(nd_ws, 'SSN'); ws_sessions{ws_count} = nd_ws.SSN;
        else; ws_sessions{ws_count} = 'unknown'; end

        if isnan(mod_N) || mod_T >= mod_N
            pref_vor = fr_T_vor; pref_spont = fr_T_spont;
            nonpref_vor = fr_N_vor; nonpref_spont = fr_N_spont;
        else
            pref_vor = fr_N_vor; pref_spont = fr_N_spont;
            nonpref_vor = fr_T_vor; nonpref_spont = fr_T_spont;
        end

        % Raw traces
        ws_pref_vor_raw(ws_count, :) = pref_vor;
        ws_pref_spont_raw(ws_count, :) = pref_spont;
        ws_nonpref_vor_raw(ws_count, :) = nonpref_vor;
        ws_nonpref_spont_raw(ws_count, :) = nonpref_spont;

        % Z-scored traces
        ref_mean_ws = mean(pref_vor(baseline_bins_mask), 'omitnan');
        ref_sd_ws = std(pref_vor(baseline_bins_mask), 'omitnan');
        if ref_sd_ws < 1; ref_sd_ws = 1; end
        ws_pref_vor(ws_count, :) = (pref_vor - ref_mean_ws) / ref_sd_ws;
        ws_pref_spont(ws_count, :) = (pref_spont - ref_mean_ws) / ref_sd_ws;
        ws_nonpref_vor(ws_count, :) = (nonpref_vor - ref_mean_ws) / ref_sd_ws;
        ws_nonpref_spont(ws_count, :) = (nonpref_spont - ref_mean_ws) / ref_sd_ws;
    end

    ws_n_sessions = length(unique(ws_sessions));
    fprintf('  %d neurons from %d sessions with >= %d saccades BOTH dirs\n', ...
        ws_count, ws_n_sessions, min_saccades);

    if ws_count < 2
        warning('Not enough neurons (%d).', ws_count);
        batch_summary = table(); all_results = precomputed_results; return
    end

    % --- Population stats ---
    fprintf('\n  --- Population PETH Statistics ---\n');
    compute_and_print_pop_stats(ws_pref_vor, ws_pref_spont, ws_nonpref_vor, ...
        ws_nonpref_spont, burst_bins_mask, baseline_bins_mask, alpha_threshold);

    % --- Z-scored figure ---
    fh_z = make_population_figure(ws_pref_vor, ws_pref_spont, ws_nonpref_vor, ...
        ws_nonpref_spont, peth_centers, ws_count, ws_n_sessions, min_saccades, ...
        'zscore', peth_window, FontSize, TickFontSize, TitleFontSize, doShow);
    print(fh_z, fullfile(out_folder, 'population_PETH_zscore.png'), '-dpng', '-r600');
    savefig(fh_z, fullfile(out_folder, 'population_PETH_zscore.fig'));
    close(fh_z);

    % --- Raw FR figure ---
    fh_raw = make_population_figure(ws_pref_vor_raw, ws_pref_spont_raw, ...
        ws_nonpref_vor_raw, ws_nonpref_spont_raw, peth_centers, ws_count, ...
        ws_n_sessions, min_saccades, 'raw', peth_window, FontSize, ...
        TickFontSize, TitleFontSize, doShow);
    print(fh_raw, fullfile(out_folder, 'population_PETH_raw.png'), '-dpng', '-r600');
    savefig(fh_raw, fullfile(out_folder, 'population_PETH_raw.fig'));
    close(fh_raw);

    % --- Best cells list ---
    print_best_cells(precomputed_results, ws_format, burst_bins_mask, ...
        baseline_bins_mask, n_best_cells, min_saccades);

    % Save
    pop_data = struct('pref_vor_z', ws_pref_vor, 'pref_spont_z', ws_pref_spont, ...
        'nonpref_vor_z', ws_nonpref_vor, 'nonpref_spont_z', ws_nonpref_spont, ...
        'pref_vor_raw', ws_pref_vor_raw, 'pref_spont_raw', ws_pref_spont_raw, ...
        'nonpref_vor_raw', ws_nonpref_vor_raw, 'nonpref_spont_raw', ws_nonpref_spont_raw, ...
        'names', {ws_names}, 'n_neurons', ws_count, 'n_sessions', ws_n_sessions, ...
        'peth_centers', peth_centers);
    save(fullfile(out_folder, 'population_PETH_data.mat'), 'pop_data', 'params');
    fprintf('\n  Saved to: %s\n', out_folder);

    batch_summary = table(); all_results = precomputed_results; params.pop_data = pop_data;
    return
end

%% ========================================================================
%  FULL RUN MODE — BUILD CELL LIST
%  ========================================================================
if isempty(cell_list)
    tfl = load(tfilelist_path, 'NPH_all_tfilelist');
    cell_list = tfl.NPH_all_tfilelist(:);
    if verbose; fprintf('Loaded %d neurons from tfilelist.\n', length(cell_list)); end
end
nNeurons = length(cell_list);
params.n_neurons = nNeurons;

% Output folder
date_str = datestr(now, 'yyyy-mm-dd');
run_num = 1;
out_folder = fullfile(save_root, sprintf('stationary_burst_analysis_%s_run%d', date_str, run_num));
while exist(out_folder, 'dir'); run_num = run_num + 1;
    out_folder = fullfile(save_root, sprintf('stationary_burst_analysis_%s_run%d', date_str, run_num)); end
mkdir(out_folder);
if doMakeFigs; fig_folder = fullfile(out_folder, 'figures'); mkdir(fig_folder); end
if verbose; fprintf('Output folder: %s\n\n', out_folder); end

%% Pre-allocate summary
sum_neuron      = cell(nNeurons,1);     sum_SSN = cell(nNeurons,1);
sum_nCCW_all    = nan(nNeurons,1);      sum_nCW_all = nan(nNeurons,1);
sum_nCCW_stat   = nan(nNeurons,1);      sum_nCW_stat = nan(nNeurons,1);
sum_nCCW_move   = nan(nNeurons,1);      sum_nCW_move = nan(nNeurons,1);
sum_nLaser_excl = nan(nNeurons,1);
sum_pref_dir    = cell(nNeurons,1);
sum_bm_vor      = nan(nNeurons,1);      sum_bm_spont = nan(nNeurons,1);
sum_p_vor       = nan(nNeurons,1);      sum_p_spont = nan(nNeurons,1);
sum_sig_vor     = false(nNeurons,1);    sum_sig_spont = false(nNeurons,1);
sum_enough      = false(nNeurons,1);
sum_stat_dur    = nan(nNeurons,1);
sum_imbalance_move = false(nNeurons,1); sum_imbalance_stat = false(nNeurons,1);
sum_status      = cell(nNeurons,1);

all_results = cell(nNeurons, 1);

%% ========================================================================
%  MAIN LOOP
%  ========================================================================
for iN = 1:nNeurons
    tfile_path = cell_list{iN};
    [sess_path, neuron_name] = fileparts(tfile_path);
    if verbose; fprintf('\n=== %d/%d: %s ===\n', iN, nNeurons, neuron_name); end

    nd = struct();
    nd.neuron_name = neuron_name;
    nd.tfile_path = tfile_path;
    nd.sess_path = sess_path;
    sum_neuron{iN} = neuron_name;

    try
        pushdir(sess_path);
        SSN = HD_GetSSN;
        nd.SSN = SSN; sum_SSN{iN} = SSN;

        % --- Load spikes ---
        [S_all, spike_cell_names, ~] = LoadSpikesQuick();
        cell_idx = find_cell_index(spike_cell_names, neuron_name, SSN);
        if isempty(cell_idx); popdir; sum_status{iN} = 'ERROR: cell not found';
            all_results{iN} = nd; continue; end
        myCell = SelectTS([], S_all, cell_idx);

        % --- Load saccade data ---
        sacc_file = strcat(SSN, '-saccades-edited.mat');
        if ~exist(sacc_file, 'file'); popdir; sum_status{iN} = 'ERROR: no saccades file';
            all_results{iN} = nd; continue; end
        sacc_data = load(sacc_file, 'temporalSaccades', 'nasalSaccades');
        temporal_all = sacc_data.temporalSaccades(~isnan(sacc_data.temporalSaccades));
        nasal_all = sacc_data.nasalSaccades(~isnan(sacc_data.nasalSaccades));

        % --- Load stationary times ---
        stat_file = strcat(SSN, '-AHV_StationaryTimes.mat');
        if ~exist(stat_file, 'file'); popdir; sum_status{iN} = 'ERROR: no StationaryTimes';
            all_results{iN} = nd; continue; end
        stat_data = load(stat_file, 'STtstart', 'STtend');
        STtstart = stat_data.STtstart; STtend = stat_data.STtend;

        % --- Separate moving vs stationary ---
        t_stat_mask = false(size(temporal_all));
        n_stat_mask = false(size(nasal_all));
        for iST = 1:length(STtstart)
            t_stat_mask = t_stat_mask | (temporal_all > STtstart(iST) & temporal_all < STtend(iST));
            n_stat_mask = n_stat_mask | (nasal_all > STtstart(iST) & nasal_all < STtend(iST));
        end
        temporal_stat = temporal_all(t_stat_mask);   nasal_stat = nasal_all(n_stat_mask);
        temporal_move = temporal_all(~t_stat_mask);  nasal_move = nasal_all(~n_stat_mask);

        % --- Laser bout exclusion (optotagged sessions only) ---
        n_laser_excluded = 0;
        is_opto_session = any(strcmp(SSN, opto_sessions));
        if is_opto_session
            try
                [~, ~, laser_on_times, ~] = SortBrainstemEventLabels3;
                if ~isnan(laser_on_times)
                    bout_windows = detect_laser_bouts(laser_on_times, laser_bout_buffer);
                    % Exclude stationary saccades within laser bout windows
                    n_before_T = length(temporal_stat);
                    n_before_N = length(nasal_stat);
                    for iBout = 1:size(bout_windows, 1)
                        temporal_stat = temporal_stat(temporal_stat < bout_windows(iBout,1) | ...
                            temporal_stat > bout_windows(iBout,2));
                        nasal_stat = nasal_stat(nasal_stat < bout_windows(iBout,1) | ...
                            nasal_stat > bout_windows(iBout,2));
                    end
                    n_laser_excluded = (n_before_T - length(temporal_stat)) + ...
                        (n_before_N - length(nasal_stat));
                    if verbose && n_laser_excluded > 0
                        fprintf('  Laser: excluded %d stationary saccades from %d bout windows\n', ...
                            n_laser_excluded, size(bout_windows, 1));
                    end
                end
            catch ME_laser
                if verbose; fprintf('  Laser detection failed: %s\n', ME_laser.message); end
            end
        end

        nT_all = length(temporal_all); nN_all = length(nasal_all);
        nT_stat = length(temporal_stat); nN_stat = length(nasal_stat);
        nT_move = length(temporal_move); nN_move = length(nasal_move);

        nd.nCCW_all = nT_all;   nd.nCW_all = nN_all;
        nd.nCCW_stat = nT_stat;  nd.nCW_stat = nN_stat;
        nd.nCCW_move = nT_move;  nd.nCW_move = nN_move;
        nd.stat_dur = sum(STtend - STtstart);
        nd.n_stat_periods = length(STtstart);
        nd.is_opto = is_opto_session;
        nd.n_laser_excluded = n_laser_excluded;

        sum_nCCW_all(iN) = nT_all;     sum_nCW_all(iN) = nN_all;
        sum_nCCW_stat(iN) = nT_stat;   sum_nCW_stat(iN) = nN_stat;
        sum_nCCW_move(iN) = nT_move;   sum_nCW_move(iN) = nN_move;
        sum_nLaser_excl(iN) = n_laser_excluded;
        sum_stat_dur(iN) = nd.stat_dur;

        % --- Saccade imbalance flags ---
        move_ratio = max(nT_move, nN_move) / max(min(nT_move, nN_move), 1);
        stat_ratio = max(nT_stat, nN_stat) / max(min(nT_stat, nN_stat), 1);
        sum_imbalance_move(iN) = move_ratio > imbalance_ratio;
        sum_imbalance_stat(iN) = stat_ratio > imbalance_ratio & (nT_stat + nN_stat) > 5;
        nd.imbalance_move = sum_imbalance_move(iN);
        nd.imbalance_stat = sum_imbalance_stat(iN);

        if verbose
            fprintf('  VOR:   CCW=%d CW=%d | Spont: CCW=%d CW=%d', nT_move, nN_move, nT_stat, nN_stat);
            if sum_imbalance_stat(iN); fprintf(' [IMBALANCED]'); end
            fprintf('\n');
        end

        % --- Compute 4 PETHs: {VOR, Spont} x {CCW, CW} ---
        cond_labels = {'move_CCW', 'move_CW', 'stat_CCW', 'stat_CW'};
        cond_sacc   = {temporal_move, nasal_move, temporal_stat, nasal_stat};
        cond_n      = {nT_move, nN_move, nT_stat, nN_stat};

        for iCond = 1:4
            if cond_n{iCond} < 1
                nd.(cond_labels{iCond}) = empty_cond_struct(nBins);
                continue
            end
            [spk_rel, evt_idx] = compute_raster_data(myCell, cond_sacc{iCond}, peth_window, peth_excessBounds);
            [m_hist, ~] = histcounts(spk_rel, peth_bins);
            fr_peth = m_hist / peth_dt / max(cond_n{iCond}, 1);
            fr_smooth = smoothdata(movmean(fr_peth, peth_smooth_win), 'gaussian', peth_smooth_win);

            [trial_burst, trial_base] = compute_per_trial_fr(myCell, cond_sacc{iCond}, burst_window, baseline_window);
            mean_burst = mean(trial_burst, 'omitnan');
            mean_base  = mean(trial_base, 'omitnan');
            bm_val = compute_burst_mod(mean_burst, mean_base);
            valid_tr = ~isnan(trial_burst) & ~isnan(trial_base);
            if sum(valid_tr) >= 3; [p_val, ~] = signrank(trial_burst(valid_tr), trial_base(valid_tr));
            else; p_val = NaN; end

            nd.(cond_labels{iCond}) = struct('fr_smooth', fr_smooth, 'burst_mod', bm_val, ...
                'p_value', p_val, 'is_sig', p_val < alpha_threshold, ...
                'mean_burst_fr', mean_burst, 'mean_baseline_fr', mean_base, ...
                'nEvents', cond_n{iCond}, 'spk_rel', spk_rel, 'evt_idx', evt_idx, ...
                'trial_burst_fr', trial_burst, 'trial_baseline_fr', trial_base);
        end

        % --- Preferred direction ---
        mod_CCW_vor = nd.move_CCW.burst_mod;
        mod_CW_vor  = nd.move_CW.burst_mod;
        if ~isempty(pref_direction_list) && isa(pref_direction_list, 'containers.Map') ...
                && pref_direction_list.isKey(neuron_name)
            pref_dir = pref_direction_list(neuron_name);
        elseif isnan(mod_CCW_vor) && isnan(mod_CW_vor); pref_dir = 'none';
        elseif isnan(mod_CW_vor) || mod_CCW_vor >= mod_CW_vor; pref_dir = 'CCW';
        else; pref_dir = 'CW';
        end
        nd.pref_direction = pref_dir;
        sum_pref_dir{iN} = pref_dir;

        if strcmp(pref_dir, 'CCW'); ps = 'CCW'; nps = 'CW';
        else; ps = 'CW'; nps = 'CCW'; end

        % Store pref/nonpref traces
        nd.pref_fr_vor    = nd.(['move_' ps]).fr_smooth;
        nd.pref_fr_spont  = nd.(['stat_' ps]).fr_smooth;
        nd.nonpref_fr_vor = nd.(['move_' nps]).fr_smooth;
        nd.nonpref_fr_spont = nd.(['stat_' nps]).fr_smooth;

        sum_bm_vor(iN)   = nd.(['move_' ps]).burst_mod;
        sum_bm_spont(iN) = nd.(['stat_' ps]).burst_mod;
        sum_p_vor(iN)    = nd.(['move_' ps]).p_value;
        sum_p_spont(iN)  = nd.(['stat_' ps]).p_value;
        sum_sig_vor(iN)  = sum_p_vor(iN) < alpha_threshold;
        sum_sig_spont(iN) = sum_p_spont(iN) < alpha_threshold;

        enough = (nT_stat >= min_saccades) && (nN_stat >= min_saccades);
        nd.enough_stationary = enough;
        sum_enough(iN) = enough;
        nd.peth_centers = peth_centers;

        % --- Per-neuron figure ---
        if doMakeFigs && enough
            fh_single = make_single_neuron_figure(nd, peth_centers, peth_window, ...
                FontSize, TickFontSize, TitleFontSize, MarkerSize_raster, LineWidth_FR, ...
                color_vor_pref, color_vor_nonpref, color_spont_pref, color_spont_nonpref, ...
                rasterColor, figWidth, figHeight, min_saccades, alpha_threshold, doShow);
            safe_name = strrep(neuron_name, '-', '_');
            if doSavePNG; print(fh_single, fullfile(fig_folder, sprintf('statBurst_%s.png', safe_name)), '-dpng', '-r600'); end
            if doSaveFIG; savefig(fh_single, fullfile(fig_folder, sprintf('statBurst_%s.fig', safe_name))); end
            close(fh_single);
        end

        sum_status{iN} = 'OK';
        popdir;
    catch ME
        sum_status{iN} = sprintf('ERROR: %s', ME.message);
        if verbose; fprintf('  ERROR: %s\n', ME.message); end
        try popdir; catch; end %#ok<CTCH>
    end
    all_results{iN} = nd;
end

%% ========================================================================
%  SUMMARY TABLE
%  ========================================================================
batch_summary = table(sum_neuron, sum_SSN, sum_nCCW_all, sum_nCW_all, ...
    sum_nCCW_move, sum_nCW_move, sum_nCCW_stat, sum_nCW_stat, sum_nLaser_excl, ...
    sum_pref_dir, sum_bm_vor, sum_bm_spont, sum_p_vor, sum_p_spont, ...
    sum_sig_vor, sum_sig_spont, sum_enough, sum_imbalance_move, sum_imbalance_stat, ...
    sum_stat_dur, sum_status, ...
    'VariableNames', {'neuron','SSN','nCCW_all','nCW_all','nCCW_move','nCW_move', ...
    'nCCW_stat','nCW_stat','nLaser_excluded','pref_dir','bm_vor','bm_spont', ...
    'p_vor','p_spont','sig_vor','sig_spont','enough_both','imbal_move','imbal_stat', ...
    'stat_dur_s','status'});

%% ========================================================================
%  PRINT SUMMARY
%  ========================================================================
n_ok = sum(strcmp(sum_status, 'OK'));
n_enough = sum(sum_enough);
n_sig_vor = sum(sum_sig_vor);
n_sig_spont = sum(sum_sig_spont & sum_enough);
n_imbal_stat = sum(sum_imbalance_stat & sum_enough);
n_opto = sum(cellfun(@(r) isfield(r,'is_opto') && r.is_opto, all_results));

fprintf('\n\n====== POPULATION SUMMARY ======\n');
fprintf('  Total neurons:                      %d\n', nNeurons);
fprintf('  Successfully processed:             %d\n', n_ok);
fprintf('  Optotagged sessions (laser excl):   %d neurons\n', n_opto);
fprintf('  With >= %d spont saccades BOTH dirs: %d neurons\n', min_saccades, n_enough);

% Threshold table
fprintf('\n  --- Threshold Table ---\n');
fprintf('  %8s  %12s  %10s\n', 'Min/dir', 'Neurons', 'Sessions');
for thr = [3 5 8 10 15 20]
    n_thr = sum(sum_nCCW_stat >= thr & sum_nCW_stat >= thr);
    s_thr = length(unique(sum_SSN(sum_nCCW_stat >= thr & sum_nCW_stat >= thr)));
    marker = ''; if thr == min_saccades; marker = ' <--'; end
    fprintf('  %8d  %12d  %10d%s\n', thr, n_thr, s_thr, marker);
end

fprintf('\n  Sig burst (VOR, pref):              %d / %d (%.0f%%)\n', n_sig_vor, n_ok, 100*n_sig_vor/max(n_ok,1));
fprintf('  Sig burst (Spont, pref):            %d / %d with enough (%.0f%%)\n', n_sig_spont, n_enough, 100*n_sig_spont/max(n_enough,1));
fprintf('  Saccade imbalance (stationary):     %d flagged\n', n_imbal_stat);

% Print imbalanced sessions
if n_imbal_stat > 0
    fprintf('\n  --- Imbalanced Stationary Saccade Sessions ---\n');
    for iN = 1:nNeurons
        if sum_imbalance_stat(iN) && sum_enough(iN)
            fprintf('    %s: CCW=%d CW=%d (ratio=%.1f)\n', sum_neuron{iN}, ...
                sum_nCCW_stat(iN), sum_nCW_stat(iN), ...
                max(sum_nCCW_stat(iN),sum_nCW_stat(iN))/max(min(sum_nCCW_stat(iN),sum_nCW_stat(iN)),1));
        end
    end
end
fprintf('================================\n');

%% ========================================================================
%  POPULATION FIGURES & STATS
%  ========================================================================
if doPopulationFig && n_enough >= 2
    % Collect data
    pop_pref_vor_z = []; pop_pref_spont_z = [];
    pop_nonpref_vor_z = []; pop_nonpref_spont_z = [];
    pop_pref_vor_raw = []; pop_pref_spont_raw = [];
    pop_nonpref_vor_raw = []; pop_nonpref_spont_raw = [];
    pop_names = {}; pop_sessions_list = {};
    pop_count = 0;

    for iN = 1:nNeurons
        if ~sum_enough(iN); continue; end
        nd_p = all_results{iN};
        if ~isfield(nd_p, 'pref_fr_vor'); continue; end
        pop_count = pop_count + 1;
        pop_names{pop_count} = nd_p.neuron_name;
        pop_sessions_list{pop_count} = nd_p.SSN;

        pop_pref_vor_raw(pop_count,:) = nd_p.pref_fr_vor;
        pop_pref_spont_raw(pop_count,:) = nd_p.pref_fr_spont;
        pop_nonpref_vor_raw(pop_count,:) = nd_p.nonpref_fr_vor;
        pop_nonpref_spont_raw(pop_count,:) = nd_p.nonpref_fr_spont;

        ref_m = mean(nd_p.pref_fr_vor(baseline_bins_mask), 'omitnan');
        ref_s = std(nd_p.pref_fr_vor(baseline_bins_mask), 'omitnan');
        if ref_s < 1; ref_s = 1; end
        pop_pref_vor_z(pop_count,:) = (nd_p.pref_fr_vor - ref_m) / ref_s;
        pop_pref_spont_z(pop_count,:) = (nd_p.pref_fr_spont - ref_m) / ref_s;
        pop_nonpref_vor_z(pop_count,:) = (nd_p.nonpref_fr_vor - ref_m) / ref_s;
        pop_nonpref_spont_z(pop_count,:) = (nd_p.nonpref_fr_spont - ref_m) / ref_s;
    end
    pop_n_sess = length(unique(pop_sessions_list));

    % Stats
    fprintf('\n  --- Population PETH Statistics (n=%d) ---\n', pop_count);
    compute_and_print_pop_stats(pop_pref_vor_z, pop_pref_spont_z, ...
        pop_nonpref_vor_z, pop_nonpref_spont_z, burst_bins_mask, baseline_bins_mask, alpha_threshold);

    % Z-scored figure
    fh_z = make_population_figure(pop_pref_vor_z, pop_pref_spont_z, ...
        pop_nonpref_vor_z, pop_nonpref_spont_z, peth_centers, pop_count, ...
        pop_n_sess, min_saccades, 'zscore', peth_window, FontSize, TickFontSize, TitleFontSize, doShow);
    if doSavePNG; print(fh_z, fullfile(out_folder, 'population_PETH_zscore.png'), '-dpng', '-r600'); end
    if doSaveFIG; savefig(fh_z, fullfile(out_folder, 'population_PETH_zscore.fig')); end
    close(fh_z);

    % Raw figure
    fh_raw = make_population_figure(pop_pref_vor_raw, pop_pref_spont_raw, ...
        pop_nonpref_vor_raw, pop_nonpref_spont_raw, peth_centers, pop_count, ...
        pop_n_sess, min_saccades, 'raw', peth_window, FontSize, TickFontSize, TitleFontSize, doShow);
    if doSavePNG; print(fh_raw, fullfile(out_folder, 'population_PETH_raw.png'), '-dpng', '-r600'); end
    if doSaveFIG; savefig(fh_raw, fullfile(out_folder, 'population_PETH_raw.fig')); end
    close(fh_raw);

    % Best cells
    print_best_cells(all_results, 'batch_v12', burst_bins_mask, baseline_bins_mask, n_best_cells, min_saccades);

    params.pop_n_neurons = pop_count;
    params.pop_n_sessions = pop_n_sess;
    params.pop_names = pop_names;
end

%% ========================================================================
%  SAVE
%  ========================================================================
save(fullfile(out_folder, 'stationary_burst_analysis.mat'), 'batch_summary', 'all_results', 'params', '-v7.3');
try writetable(batch_summary, fullfile(out_folder, 'batch_summary.csv')); catch; end
if verbose; fprintf('\nAll saved to: %s\n', out_folder); end
end

%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================

function fmt = detect_results_format(results)
    for ii = 1:length(results)
        nd = results{ii};
        if ~isstruct(nd); continue; end
        if isfield(nd, 'move_CCW'); fmt = 'batch_v12'; return; end
        if isfield(nd, 'all_CCW');  fmt = 'batch_v10'; return; end
        if isfield(nd, 'fr_T_all_smooth'); fmt = 'diagnostic'; return; end
    end
    fmt = 'unknown';
end

function [nT_s, nN_s, fr_T_vor, fr_N_vor, fr_T_spont, fr_N_spont, valid] = ...
        extract_fr_traces(nd, fmt)
    valid = false; nT_s = 0; nN_s = 0;
    fr_T_vor = []; fr_N_vor = []; fr_T_spont = []; fr_N_spont = [];
    switch fmt
        case 'batch_v12'
            if ~isfield(nd, 'nCCW_stat') || ~isfield(nd, 'move_CCW'); return; end
            nT_s = nd.nCCW_stat; nN_s = nd.nCW_stat;
            fr_T_vor = nd.move_CCW.fr_smooth; fr_N_vor = nd.move_CW.fr_smooth;
            fr_T_spont = nd.stat_CCW.fr_smooth; fr_N_spont = nd.stat_CW.fr_smooth;
            valid = true;
        case 'batch_v10'
            if ~isfield(nd, 'nCCW_stat') || ~isfield(nd, 'move_CCW'); return; end
            nT_s = nd.nCCW_stat; nN_s = nd.nCW_stat;
            fr_T_vor = nd.move_CCW.fr_smooth; fr_N_vor = nd.move_CW.fr_smooth;
            fr_T_spont = nd.stat_CCW.fr_smooth; fr_N_spont = nd.stat_CW.fr_smooth;
            valid = true;
        case 'diagnostic'
            if ~isfield(nd, 'nTemporal_stat') || ~isfield(nd, 'fr_T_all_smooth') ...
                    || ~isfield(nd, 'fr_T_stat_smooth') || ~isfield(nd, 'fr_N_stat_smooth')
                return
            end
            nT_s = nd.nTemporal_stat; nN_s = nd.nNasal_stat;
            % Diagnostic format has all, not move-only. Use all as VOR approximation.
            fr_T_vor = nd.fr_T_all_smooth; fr_N_vor = nd.fr_N_all_smooth;
            fr_T_spont = nd.fr_T_stat_smooth; fr_N_spont = nd.fr_N_stat_smooth;
            valid = true;
    end
end

function bm = compute_burst_mod_from_trace(fr_trace, bl_bins, burst_bins)
    bl_mean = mean(fr_trace(bl_bins), 'omitnan');
    burst_mean = mean(fr_trace(burst_bins), 'omitnan');
    if bl_mean > 0; bm = (burst_mean - bl_mean) / bl_mean; else; bm = NaN; end
end

function bm = compute_burst_mod(mean_burst, mean_base)
    if mean_base > 0; bm = (mean_burst - mean_base) / mean_base; else; bm = NaN; end
end

function s = empty_cond_struct(nBins)
    s = struct('fr_smooth', nan(1,nBins), 'burst_mod', NaN, 'p_value', NaN, ...
        'is_sig', false, 'mean_burst_fr', NaN, 'mean_baseline_fr', NaN, ...
        'nEvents', 0, 'spk_rel', [], 'evt_idx', [], ...
        'trial_burst_fr', [], 'trial_baseline_fr', []);
end

function bout_windows = detect_laser_bouts(laser_on, buffer)
% Detect laser bout windows from laser_on times.
% Bouts are groups of pulses with inter-pulse intervals ~2s.
% Returns Nx2 matrix of [bout_start, bout_end].
    if isnan(laser_on) || isempty(laser_on); bout_windows = []; return; end
    gaps = diff(laser_on);
    bout_breaks = find(gaps > 5);  % >5s gap = new bout
    bout_starts_idx = [1; bout_breaks + 1];
    bout_ends_idx   = [bout_breaks; length(laser_on)];
    nBouts = length(bout_starts_idx);
    bout_windows = zeros(nBouts, 2);
    for iB = 1:nBouts
        bout_windows(iB, 1) = laser_on(bout_starts_idx(iB)) - buffer;
        bout_windows(iB, 2) = laser_on(bout_ends_idx(iB)) + 1.0 + buffer;  % +1s for pulse duration
    end
end

function [spk_times_rel, event_indices] = compute_raster_data(myCell, saccade_times, window, excessBounds)
    spk_times_rel = []; event_indices = [];
    for iE = 1:length(saccade_times)
        S0 = restrict(myCell, saccade_times(iE)+window(1)-excessBounds, saccade_times(iE)+window(2)+excessBounds);
        S0 = restrict(S0, saccade_times(iE)+window(1), saccade_times(iE)+window(2));
        if ~isempty(S0.t{1})
            event_indices = [event_indices; repmat(iE, length(S0.t{1}), 1)]; %#ok<AGROW>
            spk_times_rel = [spk_times_rel; S0.t{1} - saccade_times(iE)]; %#ok<AGROW>
        end
    end
end

function [trial_burst_fr, trial_baseline_fr] = compute_per_trial_fr(myCell, sacc_times, burst_win, base_win)
    nE = length(sacc_times);
    trial_burst_fr = nan(nE,1); trial_baseline_fr = nan(nE,1);
    bd = diff(burst_win); bld = diff(base_win);
    for iE = 1:nE
        S_b = restrict(myCell, sacc_times(iE)+burst_win(1), sacc_times(iE)+burst_win(2));
        trial_burst_fr(iE) = length(S_b.t{1}) / bd;
        S_bl = restrict(myCell, sacc_times(iE)+base_win(1), sacc_times(iE)+base_win(2));
        trial_baseline_fr(iE) = length(S_bl.t{1}) / bld;
    end
end

function cell_idx = find_cell_index(spike_cell_names, neuron_name, SSN)
    cell_idx = [];
    for cc = 1:length(spike_cell_names)
        if strcmp(spike_cell_names{cc}, neuron_name); cell_idx = cc; return; end
    end
    for cc = 1:length(spike_cell_names)
        if contains(spike_cell_names{cc}, strrep(neuron_name, [SSN '-'], '')); cell_idx = cc; return; end
    end
end

function compute_and_print_pop_stats(pref_vor, pref_spont, nonpref_vor, nonpref_spont, ...
        burst_bins, base_bins, alpha)
    n = size(pref_vor, 1);
    % Per-neuron mean FR in burst window
    pv_burst = mean(pref_vor(:, burst_bins), 2, 'omitnan');
    ps_burst = mean(pref_spont(:, burst_bins), 2, 'omitnan');
    nv_burst = mean(nonpref_vor(:, burst_bins), 2, 'omitnan');
    ns_burst = mean(nonpref_spont(:, burst_bins), 2, 'omitnan');

    % Test 1: Preferred vs Non-preferred (VOR)
    [p1, ~] = signrank(pv_burst, nv_burst);
    fprintf('    Pref vs NonPref (VOR burst window):  p = %.4e  %s\n', p1, sig_str(p1, alpha));
    % Test 2: Preferred vs Non-preferred (Spontaneous)
    [p2, ~] = signrank(ps_burst, ns_burst);
    fprintf('    Pref vs NonPref (Spont burst window): p = %.4e  %s\n', p2, sig_str(p2, alpha));
    % Test 3: VOR vs Spontaneous (Preferred)
    [p3, ~] = signrank(pv_burst, ps_burst);
    fprintf('    VOR vs Spont (Preferred burst):       p = %.4e  %s\n', p3, sig_str(p3, alpha));
    % Test 4: VOR vs Spontaneous (Non-preferred)
    [p4, ~] = signrank(nv_burst, ns_burst);
    fprintf('    VOR vs Spont (NonPref burst):         p = %.4e  %s\n', p4, sig_str(p4, alpha));
end

function s = sig_str(p, alpha)
    if p < alpha/10; s = '***'; elseif p < alpha; s = '*'; else; s = 'n.s.'; end
end

function print_best_cells(all_results, fmt, burst_bins, base_bins, n_top, min_sacc)
    fprintf('\n  --- Top %d Cells for Spontaneous Saccade Bursting ---\n', n_top);
    fprintf('  %-35s  %6s  %6s  %8s  %8s\n', 'Neuron', 'CCW_sp', 'CW_sp', 'BurstMod', 'Score');
    fprintf('  %s\n', repmat('-', 1, 72));

    scores = []; names = {}; ccw_counts = []; cw_counts = []; bm_vals = [];
    for iN = 1:length(all_results)
        nd = all_results{iN};
        if ~isstruct(nd); continue; end
        [nT_s, nN_s, ~, ~, fr_T_spont, fr_N_spont, valid] = extract_fr_traces(nd, fmt);
        if ~valid || nT_s < min_sacc || nN_s < min_sacc; continue; end

        % Burst mod on preferred side
        bm_T = compute_burst_mod_from_trace(fr_T_spont, base_bins, burst_bins);
        bm_N = compute_burst_mod_from_trace(fr_N_spont, base_bins, burst_bins);
        bm_pref = max(bm_T, bm_N);
        min_count = min(nT_s, nN_s);
        score = bm_pref * min_count;

        scores(end+1) = score; names{end+1} = nd.neuron_name;
        ccw_counts(end+1) = nT_s; cw_counts(end+1) = nN_s; bm_vals(end+1) = bm_pref;
    end

    [~, sort_idx] = sort(scores, 'descend');
    for ii = 1:min(n_top, length(sort_idx))
        jj = sort_idx(ii);
        fprintf('  %-35s  %6d  %6d  %8.2f  %8.1f\n', ...
            names{jj}, ccw_counts(jj), cw_counts(jj), bm_vals(jj), scores(jj));
    end
end

function fh = make_population_figure(pref_vor, pref_spont, nonpref_vor, nonpref_spont, ...
        peth_centers, n_neurons, n_sessions, min_saccades, norm_method, peth_window, ...
        FontSize, TickFontSize, TitleFontSize, doShow)

    data_mat = {pref_vor, nonpref_vor; pref_spont, nonpref_spont};
    pop_mean = cell(2,2); pop_sem = cell(2,2);
    for iR = 1:2; for iC = 1:2
        pop_mean{iR,iC} = mean(data_mat{iR,iC}, 1, 'omitnan');
        pop_sem{iR,iC} = std(data_mat{iR,iC}, 0, 1, 'omitnan') / sqrt(n_neurons);
    end; end

    if doShow; fh = figure('Position', [50 50 1100 800]);
    else; fh = figure('Position', [50 50 1100 800], 'Visible', 'off'); end
    fh.UserData.created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fh.UserData.function = 'analyze_stationary_saccade_bursts_ver_1_2';

    row_labels = {'VOR Saccades', 'Spontaneous Saccades'};
    col_labels = {'Preferred', 'Non-Preferred'};
    line_c = {[0.20 0.50 0.80], [0.90 0.40 0.15]; [0.08 0.25 0.50], [0.60 0.22 0.05]};
    shade_c = {[0.72 0.84 0.96], [0.97 0.80 0.65]; [0.30 0.45 0.65], [0.78 0.52 0.38]};
    shade_alpha = [0.4, 0.4; 0.35, 0.35];

    all_ylim = [];
    for iR=1:2; for iC=1:2
        all_ylim = [all_ylim, pop_mean{iR,iC}+pop_sem{iR,iC}, pop_mean{iR,iC}-pop_sem{iR,iC}]; end; end
    y_lo = min(all_ylim) - 0.5; y_hi = max(all_ylim) + 0.5;

    lm=0.10; rm=0.05; tm=0.12; bm=0.10; rg=0.10; cg=0.08;
    th=1-tm-bm; rh=(th-rg)/2; tw=1-lm-rm; cw_ax=(tw-cg)/2;

    for iR=1:2; for iC=1:2
        xp = lm+(iC-1)*(cw_ax+cg); yp = 1-tm-iR*rh-(iR-1)*rg;
        ax = axes('Position', [xp yp cw_ax rh]); hold(ax,'on'); %#ok<LAXES>
        m_tr = pop_mean{iR,iC}; s_tr = pop_sem{iR,iC};
        fill(ax, [peth_centers fliplr(peth_centers)], [m_tr+s_tr fliplr(m_tr-s_tr)], ...
            shade_c{iR,iC}, 'EdgeColor','none', 'FaceAlpha', shade_alpha(iR,iC));
        plot(ax, peth_centers, m_tr, '-', 'Color', line_c{iR,iC}, 'LineWidth', 4);
        xline(ax, 0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
        xlim(ax, peth_window); ylim(ax, [y_lo y_hi]); ax.FontSize = TickFontSize;
        if iR==2; xlabel(ax, 'Time (s)', 'FontSize', FontSize); end
        if iC==1
            if strcmp(norm_method, 'zscore'); ylabel(ax, 'FR (z-score)', 'FontSize', FontSize);
            elseif strcmp(norm_method, 'raw'); ylabel(ax, 'FR (Hz)', 'FontSize', FontSize);
            else; ylabel(ax, 'FR', 'FontSize', FontSize); end
        end
        title(ax, sprintf('%s — %s', row_labels{iR}, col_labels{iC}), ...
            'FontSize', TitleFontSize, 'FontWeight', 'bold');
        box(ax, 'off');
    end; end

    if strcmp(norm_method, 'raw'); norm_label = 'raw FR';
    else; norm_label = norm_method; end
    sgtitle(sprintf('Population Saccade PETH — %s (%d neurons, %d sessions, >= %d saccades/dir)', ...
        norm_label, n_neurons, n_sessions, min_saccades), 'FontSize', TitleFontSize+2, 'FontWeight', 'bold');
end

function fh = make_single_neuron_figure(nd, peth_centers, peth_window, ...
        FontSize, TickFontSize, TitleFontSize, MarkerSize_raster, LineWidth_FR, ...
        color_vor_pref, color_vor_nonpref, color_spont_pref, color_spont_nonpref, ...
        rasterColor, figWidth, figHeight, min_saccades, alpha_thresh, doShow)

    if doShow; fh = figure('Position', [50 50 figWidth figHeight]);
    else; fh = figure('Position', [50 50 figWidth figHeight], 'Visible', 'off'); end
    fh.UserData.created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fh.UserData.function = 'analyze_stationary_saccade_bursts_ver_1_2';
    fh.UserData.neuron = nd.neuron_name;

    row_labels = {'VOR Saccades', 'Spontaneous'};
    col_labels = {'CCW (Temporal)', 'CW (Nasal)'};
    cond_keys = {'move_CCW','move_CW'; 'stat_CCW','stat_CW'};
    fr_colors = {color_vor_pref, color_vor_nonpref; color_spont_pref, color_spont_nonpref};

    all_fr_max = 0;
    for iR=1:2; for iC=1:2
        key = cond_keys{iR,iC};
        if isfield(nd,key) && ~all(isnan(nd.(key).fr_smooth))
            all_fr_max = max(all_fr_max, max(nd.(key).fr_smooth)); end
    end; end
    fr_ylim = ceil(all_fr_max/10)*10; if fr_ylim==0; fr_ylim=50; end

    lm=0.10; rm=0.08; tm=0.12; bm=0.10; rg=0.08; cg=0.06;
    th=1-tm-bm; rh=(th-rg)/2; tw=1-lm-rm; cw_ax=(tw-cg)/2;

    for iR=1:2; for iC=1:2
        xp=lm+(iC-1)*(cw_ax+cg); yp=1-tm-iR*rh-(iR-1)*rg;
        ax=axes('Position',[xp yp cw_ax rh]); hold(ax,'on'); %#ok<LAXES>
        cdata=nd.(cond_keys{iR,iC});
        if cdata.nEvents >= min_saccades && ~isempty(cdata.spk_rel)
            yyaxis(ax,'left');
            plot(ax,cdata.spk_rel,cdata.evt_idx+0.5,'.','Color',rasterColor,'MarkerSize',MarkerSize_raster);
            ylim(ax,[0 cdata.nEvents+1]); ax.YAxis(1).Color=[0 0 0];
            if iC==1; ylabel(ax,'Trial','FontSize',FontSize); else; set(ax.YAxis(1),'TickLabels',[]); end
            yyaxis(ax,'right');
            plot(ax,peth_centers,cdata.fr_smooth,'-','Color',fr_colors{iR,iC},'LineWidth',LineWidth_FR);
            ylim(ax,[0 fr_ylim]); ax.YAxis(2).Color=fr_colors{iR,iC};
            if iC==2; ylabel(ax,'FR (Hz)','FontSize',FontSize); else; set(ax.YAxis(2),'TickLabels',[]); end
        elseif cdata.nEvents>0
            text(ax,0.5,0.5,sprintf('n=%d\n(<%d)',cdata.nEvents,min_saccades),...
                'Units','normalized','HorizontalAlignment','center','FontSize',FontSize,'Color',[.6 .6 .6]);
        else
            text(ax,0.5,0.5,'n=0','Units','normalized','HorizontalAlignment','center','FontSize',FontSize,'Color',[.6 .6 .6]);
        end
        xline(ax,0,'--','Color',[.3 .3 .3],'LineWidth',1.5);
        xlim(ax,peth_window); ax.XAxis.FontSize=TickFontSize;
        if iR==2; xlabel(ax,'Time (s)','FontSize',FontSize); else; set(ax,'XTickLabel',[]); end
        sig_s=''; if cdata.nEvents>=min_saccades && cdata.p_value<alpha_thresh; sig_s=' *'; end
        title(ax,sprintf('%s (n=%d)%s',col_labels{iC},cdata.nEvents,sig_s),'FontSize',TitleFontSize,'FontWeight','bold');
        box(ax,'off');
    end; end

    ax_lab=axes('Position',[0 0 1 1],'Visible','off');
    for iR=1:2
        yc=1-tm-(iR-0.5)*rh-(iR-1)*rg;
        text(ax_lab,0.02,yc,row_labels{iR},'FontSize',FontSize,'FontWeight','bold',...
            'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90,'Units','normalized');
    end

    ps=nd.pref_direction;
    if strcmp(ps,'CCW'); bm_v=nd.move_CCW.burst_mod; bm_s_val=nd.stat_CCW.burst_mod;
    elseif strcmp(ps,'CW'); bm_v=nd.move_CW.burst_mod; bm_s_val=nd.stat_CW.burst_mod;
    else; bm_v=NaN; bm_s_val=NaN; end
    sgtitle(sprintf('%s — Pref:%s | BM: %.2f(VOR) / %.2f(Spont)', ...
        strrep(nd.neuron_name,'_','\_'), ps, bm_v, bm_s_val), 'FontSize',TitleFontSize+2,'FontWeight','bold');
end

%% Version History
% 2026-04-04  v1.2  Major refactor:
%                    - VOR (moving) vs Spontaneous (stationary) comparison
%                      instead of All vs Stationary
%                    - Laser bout exclusion for optotagged sessions
%                      (detect bouts from SortBrainstemEventLabels3,
%                      exclude saccades within bout_buffer of laser pulses)
%                    - Two population figures: z-scored + raw FR
%                    - Population stats: signrank tests for pref vs nonpref
%                      and VOR vs spontaneous in burst window
%                    - Best cells list (ranked by burst_mod * saccade_count)
%                    - Saccade imbalance flagging (ratio > imbalance_ratio)
%                    - Threshold table in command window output
%                    - Workspace mode auto-detects v1.0/v1.1/v1.2 formats
%                    - Reduced font sizes (18/14/16 from 22/16/20)
%                    - More distinct shading between VOR and spontaneous rows
% 2026-04-04  v1.1  min_saccades=8, AND criterion, pref_direction_list,
%                    population PETH, workspace mode
% 2026-04-03  v1.0  Initial version
