function [batch_summary, all_results, params] = analyze_stationary_saccade_bursts_ver_1_1(varargin)
% analyze_stationary_saccade_bursts_ver_1_1  Population analysis of
% saccade-evoked bursting during platform-stationary vs. moving periods.
%
%   For each neuron in the input list (default: 104 NPH neurons), this
%   function:
%     1. Loads spikes, saccade times, and stationary period annotations
%     2. Separates saccades into moving vs. stationary
%     3. Computes saccade-aligned FR PETHs for each condition x direction
%     4. Determines the preferred burst direction (from all-saccade data,
%        or from user-supplied pref_direction_list)
%     5. Quantifies burst modulation and tests significance via Wilcoxon
%        signed-rank test (per-trial burst window FR vs. baseline FR)
%     6. Generates a diagnostic figure per neuron (optional)
%     7. Generates a population-average PETH figure (preferred vs.
%        non-preferred, all vs. stationary saccades)
%     8. Saves all results + figures to a single output folder
%
%   INPUTS:
%       varargin    - name-value pairs (see defaults below)
%
%   KEY OPTIONS:
%       cell_list           - cell array of full t-file paths (default:
%                             loads NPH_all_tfilelist from tfilelist_path)
%       tfilelist_path      - path to NPH all tfilelist.mat
%       data_root           - root data directory (fallback scan for .t)
%       save_root           - analysis output path
%       pref_direction_list - containers.Map: neuron_name -> 'CCW'|'CW'.
%                             If empty, preferred direction is determined
%                             from burst modulation in all-saccade data.
%       min_saccades        - minimum stationary saccades required in BOTH
%                             CCW and CW directions (default: 8)
%       burst_window        - [start end] s (default: [-0.05 0.05])
%       baseline_window     - [start end] s (default: [-0.2 -0.1])
%       alpha_threshold     - significance level (default: 0.05)
%       peth_window         - full PETH window (default: [-0.2 0.2])
%       peth_dt             - PETH bin width in s (default: 0.005)
%       peth_smooth_win     - smoothdata gaussian window in bins (default: 5)
%       normalize_method    - 'zscore' or 'fold_change' (default: 'zscore')
%       doMakeFigs          - generate per-neuron figures (default: 1)
%       doPopulationFig     - generate population PETH figure (default: 1)
%       doSavePNG, doSaveFIG, doShow
%
%   OUTPUTS:
%       batch_summary   - MATLAB table with one row per neuron
%       all_results     - cell array of per-neuron result structs
%       params          - struct of parameter values used
%
%   USAGE:
%       % Full run on all 104 neurons:
%       [summary, results, par] = analyze_stationary_saccade_bursts_ver_1_1();
%
%       % Workspace mode — population PETH from existing all_results:
%       [~, ~, par] = analyze_stationary_saccade_bursts_ver_1_1(...
%           'precomputed_results', all_results, 'min_saccades', 10, 'doShow', 1);
%
%       % With custom preferred directions:
%       pref_map = containers.Map();
%       pref_map('M079-2020-10-28-TT07_2') = 'CCW';
%       [s, r, p] = analyze_stationary_saccade_bursts_ver_1_1('pref_direction_list', pref_map);
%
%   DEPENDENCIES:
%       LoadSpikesQuick, SelectTS, restrict, pushdir, popdir,
%       HD_GetSSN, process_varargin

%% ========================================================================
%  DEFAULTS
%  ========================================================================
cell_list           = {};
tfilelist_path      = 'C:\Jeff\U01\datatouse\NPH all tfilelist.mat';
data_root           = 'C:\Jeff\U01\datatouse';
save_root           = 'C:\Jeff\U01\datatouse\___analysis_output';

% Workspace mode: pass all_results from a previous run to skip data loading
% and go straight to population PETH generation
precomputed_results = {};    % cell array of per-neuron result structs

% Preferred direction override (empty = auto-detect from burst modulation)
pref_direction_list = [];   % containers.Map or empty

% Analysis windows
burst_window        = [-0.05 0.05];
baseline_window     = [-0.2 -0.1];
alpha_threshold     = 0.05;

% Saccade PETH parameters
peth_window         = [-0.2 0.2];
peth_dt             = 0.005;
peth_smooth_win     = 5;
peth_excessBounds   = 1;

% Stationary saccade filter — requires BOTH directions to meet threshold
min_saccades        = 8;

% Population figure normalization
normalize_method    = 'zscore';     % 'zscore' or 'fold_change'

% Figure parameters
FontSize            = 22;
TickFontSize        = 16;
TitleFontSize       = 20;
MarkerSize_raster   = 3;
LineWidth_FR        = 4;
color_all_CCW       = [0.30 0.69 0.29];
color_all_CW        = [0.90 0.45 0.15];
color_stat_CCW      = [0.20 0.55 0.20];
color_stat_CW       = [0.75 0.35 0.10];
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
%  PARAMS (analysis-relevant free parameters only)
%  ========================================================================
params = struct();
params.function_name        = 'analyze_stationary_saccade_bursts';
params.version              = '1.1';
params.date_run             = datestr(now, 'yyyy-mm-dd HH:MM:SS');
params.burst_window         = burst_window;
params.baseline_window      = baseline_window;
params.alpha_threshold      = alpha_threshold;
params.peth_window          = peth_window;
params.peth_dt              = peth_dt;
params.peth_smooth_win      = peth_smooth_win;
params.peth_excessBounds    = peth_excessBounds;
params.min_saccades         = min_saccades;
params.normalize_method     = normalize_method;
params.enough_stat_criterion = 'AND (both directions >= min_saccades)';
params.pref_direction_source = 'auto (burst modulation)';
if ~isempty(pref_direction_list)
    params.pref_direction_source = 'user-supplied pref_direction_list';
end

%% ========================================================================
%  WORKSPACE MODE — skip data loading, compute population PETH from
%  precomputed all_results
%  ========================================================================
if ~isempty(precomputed_results)
    if verbose; fprintf('=== WORKSPACE MODE: using precomputed_results (%d entries) ===\n', length(precomputed_results)); end

    % Create output folder
    date_str = datestr(now, 'yyyy-mm-dd');
    base_name = sprintf('stationary_burst_analysis_%s', date_str);
    run_num = 1;
    out_folder = fullfile(save_root, sprintf('%s_run%d', base_name, run_num));
    while exist(out_folder, 'dir')
        run_num = run_num + 1;
        out_folder = fullfile(save_root, sprintf('%s_run%d', base_name, run_num));
    end
    mkdir(out_folder);
    if verbose; fprintf('Output folder: %s\n\n', out_folder); end

    % Detect struct format from first fully-completed entry
    sample_nd = [];
    for ii = 1:length(precomputed_results)
        nd_test = precomputed_results{ii};
        if ~isstruct(nd_test); continue; end
        if isfield(nd_test, 'all_CCW')
            sample_nd = nd_test; break
        elseif isfield(nd_test, 'fr_T_all_smooth')
            sample_nd = nd_test; break
        end
    end
    if isempty(sample_nd); error('No fully-computed structs in precomputed_results.'); end

    if isfield(sample_nd, 'all_CCW')
        ws_format = 'batch';
    elseif isfield(sample_nd, 'fr_T_all_smooth')
        ws_format = 'diagnostic';
    else
        error('Unrecognized struct format in precomputed_results.');
    end
    if verbose; fprintf('  Detected format: %s\n', ws_format); end

    % Get peth_centers from first valid result
    ws_peth_centers = [];
    for ii = 1:length(precomputed_results)
        if isfield(precomputed_results{ii}, 'peth_centers') && ~isempty(precomputed_results{ii}.peth_centers)
            ws_peth_centers = precomputed_results{ii}.peth_centers;
            break
        end
    end
    if isempty(ws_peth_centers); error('No peth_centers found in precomputed_results.'); end
    ws_nBins = length(ws_peth_centers);
    ws_baseline_bins = ws_peth_centers >= baseline_window(1) & ws_peth_centers <= baseline_window(2);
    ws_burst_bins    = ws_peth_centers >= burst_window(1) & ws_peth_centers <= burst_window(2);

    % Collect neurons meeting threshold (both dirs >= min_saccades)
    ws_pref_all    = [];
    ws_pref_stat   = [];
    ws_nonpref_all = [];
    ws_nonpref_stat = [];
    ws_names       = {};
    ws_count       = 0;

    for iN = 1:length(precomputed_results)
        nd_ws = precomputed_results{iN};
        if ~isstruct(nd_ws); continue; end

        % Extract counts and FR traces based on format
        if strcmp(ws_format, 'batch')
            if ~isfield(nd_ws, 'nCCW_stat'); continue; end
            nT_ws = nd_ws.nCCW_stat;   nN_ws = nd_ws.nCW_stat;
            fr_T_all_ws  = nd_ws.all_CCW.fr_smooth;
            fr_N_all_ws  = nd_ws.all_CW.fr_smooth;
            fr_T_stat_ws = nd_ws.stat_CCW.fr_smooth;
            fr_N_stat_ws = nd_ws.stat_CW.fr_smooth;
            mod_T_ws = nd_ws.all_CCW.burst_mod;
            mod_N_ws = nd_ws.all_CW.burst_mod;
        else  % diagnostic
            if ~isfield(nd_ws, 'nTemporal_stat') || ~isfield(nd_ws, 'fr_T_all_smooth') ...
                    || ~isfield(nd_ws, 'fr_T_stat_smooth') || ~isfield(nd_ws, 'fr_N_stat_smooth')
                continue
            end
            nT_ws = nd_ws.nTemporal_stat;   nN_ws = nd_ws.nNasal_stat;
            fr_T_all_ws  = nd_ws.fr_T_all_smooth;
            fr_N_all_ws  = nd_ws.fr_N_all_smooth;
            fr_T_stat_ws = nd_ws.fr_T_stat_smooth;
            fr_N_stat_ws = nd_ws.fr_N_stat_smooth;
            % Compute burst modulation on the fly
            base_T_ws = mean(fr_T_all_ws(ws_baseline_bins), 'omitnan');
            base_N_ws = mean(fr_N_all_ws(ws_baseline_bins), 'omitnan');
            peak_T_ws = mean(fr_T_all_ws(ws_burst_bins), 'omitnan');
            peak_N_ws = mean(fr_N_all_ws(ws_burst_bins), 'omitnan');
            if base_T_ws > 0; mod_T_ws = (peak_T_ws - base_T_ws) / base_T_ws; else; mod_T_ws = NaN; end
            if base_N_ws > 0; mod_N_ws = (peak_N_ws - base_N_ws) / base_N_ws; else; mod_N_ws = NaN; end
        end

        % Require BOTH directions >= min_saccades
        if nT_ws < min_saccades || nN_ws < min_saccades; continue; end

        % Check valid FR data
        if all(isnan(fr_T_all_ws)) || all(isnan(fr_N_all_ws)); continue; end
        if all(isnan(fr_T_stat_ws)) || all(isnan(fr_N_stat_ws)); continue; end

        % Determine preferred direction
        if isnan(mod_T_ws) && isnan(mod_N_ws); continue; end
        if ~isempty(pref_direction_list) && isa(pref_direction_list, 'containers.Map') ...
                && pref_direction_list.isKey(nd_ws.neuron_name)
            ws_pref = pref_direction_list(nd_ws.neuron_name);
            if strcmp(ws_pref, 'CCW')
                mod_T_ws = 1; mod_N_ws = 0;  % force CCW preferred
            else
                mod_T_ws = 0; mod_N_ws = 1;  % force CW preferred
            end
        end

        if isnan(mod_N_ws) || mod_T_ws >= mod_N_ws
            fr_pref_all_ws  = fr_T_all_ws;   fr_pref_stat_ws  = fr_T_stat_ws;
            fr_nonpref_all_ws = fr_N_all_ws;  fr_nonpref_stat_ws = fr_N_stat_ws;
        else
            fr_pref_all_ws  = fr_N_all_ws;   fr_pref_stat_ws  = fr_N_stat_ws;
            fr_nonpref_all_ws = fr_T_all_ws;  fr_nonpref_stat_ws = fr_T_stat_ws;
        end

        ws_count = ws_count + 1;
        ws_names{ws_count} = nd_ws.neuron_name;

        % Normalize
        ref_mean_ws = mean(fr_pref_all_ws(ws_baseline_bins), 'omitnan');
        ref_sd_ws   = std(fr_pref_all_ws(ws_baseline_bins), 'omitnan');

        if strcmp(normalize_method, 'zscore')
            if ref_sd_ws < 1; ref_sd_ws = 1; end
            ws_pref_all(ws_count, :)     = (fr_pref_all_ws - ref_mean_ws) / ref_sd_ws;
            ws_pref_stat(ws_count, :)    = (fr_pref_stat_ws - ref_mean_ws) / ref_sd_ws;
            ws_nonpref_all(ws_count, :)  = (fr_nonpref_all_ws - ref_mean_ws) / ref_sd_ws;
            ws_nonpref_stat(ws_count, :) = (fr_nonpref_stat_ws - ref_mean_ws) / ref_sd_ws;
        else
            if ref_mean_ws < 1; ref_mean_ws = 1; end
            ws_pref_all(ws_count, :)     = fr_pref_all_ws / ref_mean_ws;
            ws_pref_stat(ws_count, :)    = fr_pref_stat_ws / ref_mean_ws;
            ws_nonpref_all(ws_count, :)  = fr_nonpref_all_ws / ref_mean_ws;
            ws_nonpref_stat(ws_count, :) = fr_nonpref_stat_ws / ref_mean_ws;
        end
    end

    fprintf('  Population PETH: %d neurons with >= %d saccades BOTH dirs\n', ws_count, min_saccades);

    % Count unique sessions
    ws_sessions = {};
    for ii = 1:ws_count
        nd_sess = [];
        for jj = 1:length(precomputed_results)
            if strcmp(precomputed_results{jj}.neuron_name, ws_names{ii})
                if isfield(precomputed_results{jj}, 'SSN')
                    nd_sess = precomputed_results{jj}.SSN;
                else
                    % Extract session from neuron name (remove TT part)
                    tok = regexp(ws_names{ii}, '(M\d+-\d+-\d+-\d+(?:-\d+)?)', 'tokens');
                    if ~isempty(tok); nd_sess = tok{1}{1}; end
                end
                break
            end
        end
        if ~isempty(nd_sess); ws_sessions{end+1} = nd_sess; end %#ok<AGROW>
    end
    ws_n_sessions = length(unique(ws_sessions));
    fprintf('  From %d unique sessions\n', ws_n_sessions);

    if ws_count < 2
        warning('Not enough neurons (%d). Try lowering min_saccades.', ws_count);
        batch_summary = table();
        all_results = precomputed_results;
        return
    end

    % Generate population figure
    fh_ws = make_population_figure(ws_pref_all, ws_pref_stat, ws_nonpref_all, ...
        ws_nonpref_stat, ws_peth_centers, ws_count, ws_n_sessions, min_saccades, ...
        normalize_method, peth_window, FontSize, TickFontSize, TitleFontSize, doShow);

    if doSavePNG
        print(fh_ws, fullfile(out_folder, 'population_PETH.png'), '-dpng', '-r600');
    end
    if doSaveFIG
        savefig(fh_ws, fullfile(out_folder, 'population_PETH.fig'));
    end

    % Save data
    pop_ws_data = struct();
    pop_ws_data.pref_all     = ws_pref_all;
    pop_ws_data.pref_stat    = ws_pref_stat;
    pop_ws_data.nonpref_all  = ws_nonpref_all;
    pop_ws_data.nonpref_stat = ws_nonpref_stat;
    pop_ws_data.neuron_names = ws_names;
    pop_ws_data.n_neurons    = ws_count;
    pop_ws_data.peth_centers = ws_peth_centers;
    pop_ws_data.min_saccades = min_saccades;
    pop_ws_data.normalize_method = normalize_method;
    save(fullfile(out_folder, 'population_PETH_data.mat'), 'pop_ws_data', 'params');
    fprintf('  Saved to: %s\n', out_folder);

    batch_summary = table();
    all_results = precomputed_results;
    params.pop_data = pop_ws_data;
    return
end

%% ========================================================================
%  BUILD CELL LIST
%  ========================================================================
if isempty(cell_list)
    if exist(tfilelist_path, 'file')
        tfl = load(tfilelist_path, 'NPH_all_tfilelist');
        cell_list = tfl.NPH_all_tfilelist(:);
        if verbose; fprintf('Loaded %d neurons from tfilelist.\n', length(cell_list)); end
    else
        warning('No tfilelist found. Scanning %s for .t files...', data_root);
        t_files = dir(fullfile(data_root, '**', '*.t'));
        cell_list = cell(length(t_files), 1);
        for ii = 1:length(t_files)
            cell_list{ii} = fullfile(t_files(ii).folder, t_files(ii).name);
        end
        if verbose; fprintf('Found %d .t files in %s\n', length(cell_list), data_root); end
    end
end

nNeurons = length(cell_list);
params.n_neurons = nNeurons;

%% ========================================================================
%  CREATE OUTPUT FOLDER
%  ========================================================================
date_str = datestr(now, 'yyyy-mm-dd');
base_name = sprintf('stationary_burst_analysis_%s', date_str);
run_num = 1;
out_folder = fullfile(save_root, sprintf('%s_run%d', base_name, run_num));
while exist(out_folder, 'dir')
    run_num = run_num + 1;
    out_folder = fullfile(save_root, sprintf('%s_run%d', base_name, run_num));
end
mkdir(out_folder);
if doMakeFigs
    fig_folder = fullfile(out_folder, 'figures');
    mkdir(fig_folder);
end
if verbose; fprintf('Output folder: %s\n\n', out_folder); end

%% ========================================================================
%  PRE-ALLOCATE SUMMARY ARRAYS
%  ========================================================================
sum_neuron_name     = cell(nNeurons, 1);
sum_SSN             = cell(nNeurons, 1);
sum_nCCW_all        = nan(nNeurons, 1);
sum_nCW_all         = nan(nNeurons, 1);
sum_nCCW_stat       = nan(nNeurons, 1);
sum_nCW_stat        = nan(nNeurons, 1);
sum_nCCW_move       = nan(nNeurons, 1);
sum_nCW_move        = nan(nNeurons, 1);
sum_pref_dir        = cell(nNeurons, 1);
sum_burst_mod_all   = nan(nNeurons, 1);
sum_burst_mod_move  = nan(nNeurons, 1);
sum_burst_mod_stat  = nan(nNeurons, 1);
sum_p_all           = nan(nNeurons, 1);
sum_p_move          = nan(nNeurons, 1);
sum_p_stat          = nan(nNeurons, 1);
sum_sig_all         = false(nNeurons, 1);
sum_sig_move        = false(nNeurons, 1);
sum_sig_stat        = false(nNeurons, 1);
sum_stat_dur        = nan(nNeurons, 1);
sum_n_stat_periods  = nan(nNeurons, 1);
sum_enough_stat     = false(nNeurons, 1);
sum_status          = cell(nNeurons, 1);

all_results = cell(nNeurons, 1);

%% ========================================================================
%  MAIN LOOP
%  ========================================================================
peth_bins    = linspace(peth_window(1), peth_window(2), diff(peth_window)/peth_dt + 1);
peth_centers = peth_bins(1:end-1) + diff(peth_bins) / 2;
nBins        = length(peth_centers);

for iN = 1:nNeurons
    tfile_path = cell_list{iN};
    [sess_path, neuron_name] = fileparts(tfile_path);

    if verbose
        fprintf('\n=== %d/%d: %s ===\n', iN, nNeurons, neuron_name);
    end

    nd = struct();
    nd.neuron_name  = neuron_name;
    nd.tfile_path   = tfile_path;
    nd.sess_path    = sess_path;
    sum_neuron_name{iN} = neuron_name;

    try
        %% ---- Navigate to session ----
        pushdir(sess_path);
        SSN = HD_GetSSN;
        nd.SSN = SSN;
        sum_SSN{iN} = SSN;

        %% ---- Load spikes ----
        [S_all, spike_cell_names, ~] = LoadSpikesQuick();
        cell_idx = find_cell_index(spike_cell_names, neuron_name, SSN);
        if isempty(cell_idx)
            popdir;
            sum_status{iN} = 'ERROR: cell not found';
            if verbose; fprintf('  SKIP: could not match cell name\n'); end
            all_results{iN} = nd;
            continue
        end
        myCell = SelectTS([], S_all, cell_idx);

        %% ---- Load saccade data ----
        sacc_file = strcat(SSN, '-saccades-edited.mat');
        if ~exist(sacc_file, 'file')
            popdir;
            sum_status{iN} = 'ERROR: no saccades-edited.mat';
            all_results{iN} = nd; continue
        end
        sacc_data = load(sacc_file, 'temporalSaccades', 'nasalSaccades');
        temporal_all = sacc_data.temporalSaccades;
        nasal_all    = sacc_data.nasalSaccades;
        temporal_all = temporal_all(~isnan(temporal_all));
        nasal_all    = nasal_all(~isnan(nasal_all));

        %% ---- Load stationary times ----
        stat_file = strcat(SSN, '-AHV_StationaryTimes.mat');
        if ~exist(stat_file, 'file')
            popdir;
            sum_status{iN} = 'ERROR: no StationaryTimes.mat';
            all_results{iN} = nd; continue
        end
        stat_data = load(stat_file, 'STtstart', 'STtend');
        STtstart = stat_data.STtstart;
        STtend   = stat_data.STtend;

        %% ---- Separate stationary vs moving saccades ----
        temporal_stat_mask = false(size(temporal_all));
        nasal_stat_mask    = false(size(nasal_all));
        for iST = 1:length(STtstart)
            temporal_stat_mask = temporal_stat_mask | ...
                (temporal_all > STtstart(iST) & temporal_all < STtend(iST));
            nasal_stat_mask = nasal_stat_mask | ...
                (nasal_all > STtstart(iST) & nasal_all < STtend(iST));
        end

        temporal_stat = temporal_all(temporal_stat_mask);
        nasal_stat    = nasal_all(nasal_stat_mask);
        temporal_move = temporal_all(~temporal_stat_mask);
        nasal_move    = nasal_all(~nasal_stat_mask);

        nT_all  = length(temporal_all);     nN_all  = length(nasal_all);
        nT_stat = length(temporal_stat);    nN_stat = length(nasal_stat);
        nT_move = length(temporal_move);    nN_move = length(nasal_move);
        stat_dur_total = sum(STtend - STtstart);

        nd.nCCW_all = nT_all;  nd.nCW_all = nN_all;
        nd.nCCW_stat = nT_stat; nd.nCW_stat = nN_stat;
        nd.nCCW_move = nT_move; nd.nCW_move = nN_move;
        nd.stat_dur = stat_dur_total;
        nd.n_stat_periods = length(STtstart);

        sum_nCCW_all(iN) = nT_all;   sum_nCW_all(iN) = nN_all;
        sum_nCCW_stat(iN) = nT_stat;  sum_nCW_stat(iN) = nN_stat;
        sum_nCCW_move(iN) = nT_move;  sum_nCW_move(iN) = nN_move;
        sum_stat_dur(iN) = stat_dur_total;
        sum_n_stat_periods(iN) = length(STtstart);

        if verbose
            fprintf('  CCW: %d all | %d move | %d stat\n', nT_all, nT_move, nT_stat);
            fprintf('  CW:  %d all | %d move | %d stat\n', nN_all, nN_move, nN_stat);
        end

        %% ---- Compute PETHs and burst metrics for 6 conditions ----
        conditions = struct( ...
            'label',     {'all_CCW','all_CW','move_CCW','move_CW','stat_CCW','stat_CW'}, ...
            'sacc_times', {temporal_all, nasal_all, temporal_move, nasal_move, temporal_stat, nasal_stat}, ...
            'nEvents',   {nT_all, nN_all, nT_move, nN_move, nT_stat, nN_stat} ...
        );

        for iCond = 1:6
            cond = conditions(iCond);
            if cond.nEvents < 1
                nd.(cond.label).fr_smooth = nan(1, nBins);
                nd.(cond.label).burst_mod = NaN;
                nd.(cond.label).p_value = NaN;
                nd.(cond.label).is_sig = false;
                nd.(cond.label).mean_burst_fr = NaN;
                nd.(cond.label).mean_baseline_fr = NaN;
                nd.(cond.label).nEvents = 0;
                nd.(cond.label).spk_rel = [];
                nd.(cond.label).evt_idx = [];
                nd.(cond.label).trial_burst_fr = [];
                nd.(cond.label).trial_baseline_fr = [];
                continue
            end

            [spk_rel, evt_idx] = compute_raster_data(myCell, cond.sacc_times, ...
                peth_window, peth_excessBounds);
            [m_hist, ~] = histcounts(spk_rel, peth_bins);
            fr_peth = m_hist / peth_dt / max(cond.nEvents, 1);
            fr_smooth = smoothdata(movmean(fr_peth, peth_smooth_win), 'gaussian', peth_smooth_win);

            [trial_burst_fr, trial_baseline_fr] = compute_per_trial_fr( ...
                myCell, cond.sacc_times, burst_window, baseline_window);

            mean_burst_fr    = mean(trial_burst_fr, 'omitnan');
            mean_baseline_fr = mean(trial_baseline_fr, 'omitnan');
            if mean_baseline_fr > 0
                burst_mod_val = (mean_burst_fr - mean_baseline_fr) / mean_baseline_fr;
            else
                burst_mod_val = NaN;
            end

            valid_trials = ~isnan(trial_burst_fr) & ~isnan(trial_baseline_fr);
            if sum(valid_trials) >= 3
                [p_val, ~] = signrank(trial_burst_fr(valid_trials), ...
                    trial_baseline_fr(valid_trials));
            else
                p_val = NaN;
            end

            nd.(cond.label).fr_smooth       = fr_smooth;
            nd.(cond.label).burst_mod       = burst_mod_val;
            nd.(cond.label).p_value         = p_val;
            nd.(cond.label).is_sig          = p_val < alpha_threshold;
            nd.(cond.label).mean_burst_fr   = mean_burst_fr;
            nd.(cond.label).mean_baseline_fr = mean_baseline_fr;
            nd.(cond.label).nEvents         = cond.nEvents;
            nd.(cond.label).spk_rel         = spk_rel;
            nd.(cond.label).evt_idx         = evt_idx;
            nd.(cond.label).trial_burst_fr  = trial_burst_fr;
            nd.(cond.label).trial_baseline_fr = trial_baseline_fr;
        end

        %% ---- Determine preferred direction ----
        if ~isempty(pref_direction_list) && pref_direction_list.isKey(neuron_name)
            pref_dir = pref_direction_list(neuron_name);
            nd.pref_source = 'user_supplied';
        else
            % Auto-detect from all-saccade burst modulation
            mod_CCW = nd.all_CCW.burst_mod;
            mod_CW  = nd.all_CW.burst_mod;
            if isnan(mod_CCW) && isnan(mod_CW)
                pref_dir = 'none';
            elseif isnan(mod_CW) || mod_CCW >= mod_CW
                pref_dir = 'CCW';
            else
                pref_dir = 'CW';
            end
            nd.pref_source = 'auto_burst_mod';
        end
        nd.pref_direction = pref_dir;
        sum_pref_dir{iN} = pref_dir;

        % Map preferred/non-preferred to CCW/CW condition labels
        if strcmp(pref_dir, 'CCW')
            pref_suffix = 'CCW'; nonpref_suffix = 'CW';
        elseif strcmp(pref_dir, 'CW')
            pref_suffix = 'CW'; nonpref_suffix = 'CCW';
        else
            pref_suffix = 'CCW'; nonpref_suffix = 'CW'; % fallback
        end

        % Store preferred/non-preferred FR traces for population averaging
        nd.pref_fr_all  = nd.(['all_' pref_suffix]).fr_smooth;
        nd.pref_fr_move = nd.(['move_' pref_suffix]).fr_smooth;
        nd.pref_fr_stat = nd.(['stat_' pref_suffix]).fr_smooth;
        nd.nonpref_fr_all  = nd.(['all_' nonpref_suffix]).fr_smooth;
        nd.nonpref_fr_move = nd.(['move_' nonpref_suffix]).fr_smooth;
        nd.nonpref_fr_stat = nd.(['stat_' nonpref_suffix]).fr_smooth;

        % Summary metrics (preferred direction)
        sum_burst_mod_all(iN)  = nd.(['all_' pref_suffix]).burst_mod;
        sum_burst_mod_move(iN) = nd.(['move_' pref_suffix]).burst_mod;
        sum_burst_mod_stat(iN) = nd.(['stat_' pref_suffix]).burst_mod;
        sum_p_all(iN)          = nd.(['all_' pref_suffix]).p_value;
        sum_p_move(iN)         = nd.(['move_' pref_suffix]).p_value;
        sum_p_stat(iN)         = nd.(['stat_' pref_suffix]).p_value;

        sum_sig_all(iN)  = sum_p_all(iN) < alpha_threshold;
        sum_sig_move(iN) = sum_p_move(iN) < alpha_threshold;
        sum_sig_stat(iN) = sum_p_stat(iN) < alpha_threshold;

        % Enough stationary = BOTH directions >= min_saccades
        enough_stat = (nT_stat >= min_saccades) && (nN_stat >= min_saccades);
        nd.enough_stationary = enough_stat;
        sum_enough_stat(iN)  = enough_stat;

        nd.peth_centers = peth_centers;

        if verbose
            fprintf('  Pref: %s (%s) | BurstMod (all/move/stat): %.2f / %.2f / %.2f\n', ...
                pref_dir, nd.pref_source, sum_burst_mod_all(iN), ...
                sum_burst_mod_move(iN), sum_burst_mod_stat(iN));
            fprintf('  p-vals (all/move/stat): %.4f / %.4f / %.4f | enough_stat: %d\n', ...
                sum_p_all(iN), sum_p_move(iN), sum_p_stat(iN), enough_stat);
        end

        %% ---- Per-neuron diagnostic figure ----
        if doMakeFigs && enough_stat
            fh = make_diagnostic_figure(nd, peth_centers, peth_window, ...
                FontSize, TickFontSize, TitleFontSize, MarkerSize_raster, ...
                LineWidth_FR, color_all_CCW, color_all_CW, color_stat_CCW, ...
                color_stat_CW, rasterColor, figWidth, figHeight, ...
                min_saccades, alpha_threshold, doShow);

            safe_name = strrep(neuron_name, '-', '_');
            if doSavePNG
                print(fh, fullfile(fig_folder, sprintf('statBurst_%s.png', safe_name)), '-dpng', '-r600');
            end
            if doSaveFIG
                savefig(fh, fullfile(fig_folder, sprintf('statBurst_%s.fig', safe_name)));
            end
            close(fh);
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
%  BUILD SUMMARY TABLE
%  ========================================================================
batch_summary = table( ...
    sum_neuron_name, sum_SSN, ...
    sum_nCCW_all, sum_nCW_all, sum_nCCW_stat, sum_nCW_stat, ...
    sum_nCCW_move, sum_nCW_move, ...
    sum_pref_dir, sum_burst_mod_all, sum_burst_mod_move, sum_burst_mod_stat, ...
    sum_p_all, sum_p_move, sum_p_stat, ...
    sum_sig_all, sum_sig_move, sum_sig_stat, ...
    sum_enough_stat, sum_stat_dur, sum_n_stat_periods, sum_status, ...
    'VariableNames', { ...
    'neuron_name', 'SSN', ...
    'nCCW_all', 'nCW_all', 'nCCW_stat', 'nCW_stat', ...
    'nCCW_move', 'nCW_move', ...
    'pref_direction', 'burst_mod_all', 'burst_mod_move', 'burst_mod_stat', ...
    'p_all', 'p_move', 'p_stat', ...
    'sig_all', 'sig_move', 'sig_stat', ...
    'enough_stat_both', 'stationary_dur_s', 'n_stat_periods', 'status'});

%% ========================================================================
%  PRINT SUMMARY
%  ========================================================================
n_ok = sum(strcmp(sum_status, 'OK'));
n_enough_stat = sum(sum_enough_stat);
n_sig_all  = sum(sum_sig_all);
n_sig_move = sum(sum_sig_move);
n_sig_stat = sum(sum_sig_stat & sum_enough_stat);

fprintf('\n\n====== POPULATION SUMMARY ======\n');
fprintf('  Total neurons:                  %d\n', nNeurons);
fprintf('  Successfully processed:         %d\n', n_ok);
fprintf('  With >= %d stat saccades BOTH dirs: %d\n', min_saccades, n_enough_stat);
fprintf('  Sig burst (all sacc, pref dir):   %d / %d (%.0f%%)\n', ...
    n_sig_all, n_ok, 100*n_sig_all/max(n_ok,1));
fprintf('  Sig burst (moving, pref dir):     %d / %d (%.0f%%)\n', ...
    n_sig_move, n_ok, 100*n_sig_move/max(n_ok,1));
fprintf('  Sig burst (stationary, pref dir): %d / %d with enough (%.0f%%)\n', ...
    n_sig_stat, n_enough_stat, 100*n_sig_stat/max(n_enough_stat,1));
fprintf('================================\n');

%% ========================================================================
%  POPULATION PETH FIGURE
%  ========================================================================
if doPopulationFig && n_enough_stat >= 2
    if verbose; fprintf('\nGenerating population PETH figure (%d neurons)...\n', n_enough_stat); end

    % Collect normalized PETHs from neurons with enough stationary saccades
    pop_pref_all    = nan(n_enough_stat, nBins);
    pop_pref_stat   = nan(n_enough_stat, nBins);
    pop_nonpref_all = nan(n_enough_stat, nBins);
    pop_nonpref_stat = nan(n_enough_stat, nBins);
    pop_neuron_names = cell(n_enough_stat, 1);

    baseline_bins_pop = peth_centers >= baseline_window(1) & peth_centers <= baseline_window(2);
    pop_idx = 0;
    for iN = 1:nNeurons
        if ~sum_enough_stat(iN); continue; end
        nd_pop = all_results{iN};
        if ~isfield(nd_pop, 'pref_fr_all'); continue; end

        pop_idx = pop_idx + 1;
        pop_neuron_names{pop_idx} = nd_pop.neuron_name;

        ref_mean_pop = mean(nd_pop.pref_fr_all(baseline_bins_pop), 'omitnan');
        ref_sd_pop   = std(nd_pop.pref_fr_all(baseline_bins_pop), 'omitnan');

        if strcmp(normalize_method, 'zscore')
            if ref_sd_pop < 1; ref_sd_pop = 1; end
            pop_pref_all(pop_idx, :)     = (nd_pop.pref_fr_all - ref_mean_pop) / ref_sd_pop;
            pop_pref_stat(pop_idx, :)    = (nd_pop.pref_fr_stat - ref_mean_pop) / ref_sd_pop;
            pop_nonpref_all(pop_idx, :)  = (nd_pop.nonpref_fr_all - ref_mean_pop) / ref_sd_pop;
            pop_nonpref_stat(pop_idx, :) = (nd_pop.nonpref_fr_stat - ref_mean_pop) / ref_sd_pop;
        else
            if ref_mean_pop < 1; ref_mean_pop = 1; end
            pop_pref_all(pop_idx, :)     = nd_pop.pref_fr_all / ref_mean_pop;
            pop_pref_stat(pop_idx, :)    = nd_pop.pref_fr_stat / ref_mean_pop;
            pop_nonpref_all(pop_idx, :)  = nd_pop.nonpref_fr_all / ref_mean_pop;
            pop_nonpref_stat(pop_idx, :) = nd_pop.nonpref_fr_stat / ref_mean_pop;
        end
    end

    pop_pref_all     = pop_pref_all(1:pop_idx, :);
    pop_pref_stat    = pop_pref_stat(1:pop_idx, :);
    pop_nonpref_all  = pop_nonpref_all(1:pop_idx, :);
    pop_nonpref_stat = pop_nonpref_stat(1:pop_idx, :);

    % Count unique sessions among included neurons
    pop_sessions = cell(pop_idx, 1);
    for kk = 1:pop_idx
        for jj = 1:nNeurons
            if strcmp(all_results{jj}.neuron_name, pop_neuron_names{kk})
                pop_sessions{kk} = all_results{jj}.SSN;
                break
            end
        end
    end
    pop_n_sessions = length(unique(pop_sessions(~cellfun('isempty', pop_sessions))));

    fh_pop = make_population_figure(pop_pref_all, pop_pref_stat, ...
        pop_nonpref_all, pop_nonpref_stat, peth_centers, pop_idx, ...
        pop_n_sessions, min_saccades, normalize_method, peth_window, ...
        FontSize, TickFontSize, TitleFontSize, doShow);

    if doSavePNG
        print(fh_pop, fullfile(out_folder, 'population_PETH.png'), '-dpng', '-r600');
    end
    if doSaveFIG
        savefig(fh_pop, fullfile(out_folder, 'population_PETH.fig'));
    end
    close(fh_pop);
    if verbose; fprintf('  Population PETH figure saved.\n'); end

    params.pop_neuron_names   = pop_neuron_names(1:pop_idx);
    params.pop_n_neurons      = pop_idx;
    params.pop_pref_all       = pop_pref_all;
    params.pop_pref_stat      = pop_pref_stat;
    params.pop_nonpref_all    = pop_nonpref_all;
    params.pop_nonpref_stat   = pop_nonpref_stat;

elseif doPopulationFig
    fprintf('\n  Not enough neurons with sufficient stationary saccades for population figure.\n');
end

%% ========================================================================
%  SAVE
%  ========================================================================
save(fullfile(out_folder, 'stationary_burst_analysis.mat'), ...
    'batch_summary', 'all_results', 'params', '-v7.3');
if verbose; fprintf('\nResults saved to: %s\n', out_folder); end

try
    writetable(batch_summary, fullfile(out_folder, 'batch_summary.csv'));
    if verbose; fprintf('Summary CSV saved.\n'); end
catch
    warning('Could not write CSV summary.');
end

end

%% ========================================================================
%  LOCAL HELPER: compute raster data
%  ========================================================================
function [spk_times_rel, event_indices] = compute_raster_data(myCell, saccade_times, window, excessBounds)
    spk_times_rel   = [];
    event_indices    = [];
    nEvents         = length(saccade_times);
    for iE = 1:nEvents
        S0 = restrict(myCell, saccade_times(iE) + window(1) - excessBounds, ...
                               saccade_times(iE) + window(2) + excessBounds);
        S0 = restrict(S0, saccade_times(iE) + window(1), ...
                          saccade_times(iE) + window(2));
        if ~isempty(S0.t{1})
            nSpk = length(S0.t{1});
            event_indices = [event_indices; repmat(iE, nSpk, 1)]; %#ok<AGROW>
            spk_times_rel = [spk_times_rel; S0.t{1} - saccade_times(iE)]; %#ok<AGROW>
        end
    end
end

%% ========================================================================
%  LOCAL HELPER: compute per-trial FR in burst and baseline windows
%  ========================================================================
function [trial_burst_fr, trial_baseline_fr] = compute_per_trial_fr( ...
        myCell, saccade_times, burst_win, baseline_win)
    nEvents = length(saccade_times);
    trial_burst_fr    = nan(nEvents, 1);
    trial_baseline_fr = nan(nEvents, 1);
    burst_dur    = diff(burst_win);
    baseline_dur = diff(baseline_win);
    for iE = 1:nEvents
        t_sacc = saccade_times(iE);
        S_burst = restrict(myCell, t_sacc + burst_win(1), t_sacc + burst_win(2));
        trial_burst_fr(iE) = length(S_burst.t{1}) / burst_dur;
        S_base = restrict(myCell, t_sacc + baseline_win(1), t_sacc + baseline_win(2));
        trial_baseline_fr(iE) = length(S_base.t{1}) / baseline_dur;
    end
end

%% ========================================================================
%  LOCAL HELPER: find cell index
%  ========================================================================
function cell_idx = find_cell_index(spike_cell_names, neuron_name, SSN)
    cell_idx = [];
    for cc = 1:length(spike_cell_names)
        if strcmp(spike_cell_names{cc}, neuron_name)
            cell_idx = cc; return
        end
    end
    for cc = 1:length(spike_cell_names)
        if contains(spike_cell_names{cc}, strrep(neuron_name, [SSN '-'], ''))
            cell_idx = cc; return
        end
    end
end

%% ========================================================================
%  LOCAL HELPER: diagnostic figure (same as v1.0)
%  ========================================================================
function fh = make_diagnostic_figure(nd, peth_centers, peth_window, ...
        FontSize, TickFontSize, TitleFontSize, MarkerSize_raster, ...
        LineWidth_FR, color_all_CCW, color_all_CW, color_stat_CCW, ...
        color_stat_CW, rasterColor, figWidth, figHeight, ...
        min_saccades, alpha_thresh, doShow)

    if doShow
        fh = figure('Position', [50 50 figWidth figHeight]);
    else
        fh = figure('Position', [50 50 figWidth figHeight], 'Visible', 'off');
    end
    fh.UserData.created  = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fh.UserData.function = 'analyze_stationary_saccade_bursts_ver_1_1';
    fh.UserData.neuron   = nd.neuron_name;

    row_labels = {'All Saccades', 'Stationary Only'};
    col_labels = {'CCW (Temporal)', 'CW (Nasal)'};
    cond_keys = {'all_CCW', 'all_CW'; 'stat_CCW', 'stat_CW'};
    fr_colors = {color_all_CCW, color_all_CW; color_stat_CCW, color_stat_CW};

    all_fr_max = 0;
    for iR = 1:2
        for iC = 1:2
            key = cond_keys{iR, iC};
            if isfield(nd, key) && ~all(isnan(nd.(key).fr_smooth))
                all_fr_max = max(all_fr_max, max(nd.(key).fr_smooth));
            end
        end
    end
    fr_ylim_upper = ceil(all_fr_max / 10) * 10;
    if fr_ylim_upper == 0; fr_ylim_upper = 50; end

    left_m = 0.10; right_m = 0.08; top_m = 0.12; bottom_m = 0.10;
    row_gap = 0.08; col_gap = 0.06;
    total_h = 1 - top_m - bottom_m;
    row_h = (total_h - row_gap) / 2;
    total_w = 1 - left_m - right_m;
    col_w = (total_w - col_gap) / 2;

    for iR = 1:2
        for iC = 1:2
            x_pos = left_m + (iC-1)*(col_w + col_gap);
            y_pos = 1 - top_m - iR*row_h - (iR-1)*row_gap;
            ax = axes('Position', [x_pos, y_pos, col_w, row_h]); %#ok<LAXES>
            hold(ax, 'on');

            key = cond_keys{iR, iC};
            cdata = nd.(key);

            if cdata.nEvents >= min_saccades && ~isempty(cdata.spk_rel)
                yyaxis(ax, 'left');
                plot(ax, cdata.spk_rel, cdata.evt_idx + 0.5, '.', ...
                    'Color', rasterColor, 'MarkerSize', MarkerSize_raster);
                ylim(ax, [0 cdata.nEvents + 1]);
                ax.YAxis(1).Color = [0 0 0];
                if iC == 1; ylabel(ax, 'Trial', 'FontSize', FontSize);
                else; set(ax.YAxis(1), 'TickLabels', []); end

                yyaxis(ax, 'right');
                plot(ax, peth_centers, cdata.fr_smooth, '-', ...
                    'Color', fr_colors{iR, iC}, 'LineWidth', LineWidth_FR);
                ylim(ax, [0 fr_ylim_upper]);
                ax.YAxis(2).Color = fr_colors{iR, iC};
                if iC == 2; ylabel(ax, 'FR (Hz)', 'FontSize', FontSize);
                else; set(ax.YAxis(2), 'TickLabels', []); end
            elseif cdata.nEvents > 0 && cdata.nEvents < min_saccades
                text(ax, 0.5, 0.5, sprintf('n = %d\n(< %d)', cdata.nEvents, min_saccades), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                    'FontSize', FontSize, 'Color', [0.6 0.6 0.6]);
            else
                text(ax, 0.5, 0.5, 'n = 0', 'Units', 'normalized', ...
                    'HorizontalAlignment', 'center', 'FontSize', FontSize, ...
                    'Color', [0.6 0.6 0.6]);
            end

            xline(ax, 0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
            xlim(ax, peth_window);
            ax.XAxis.FontSize = TickFontSize;
            if iR == 2; xlabel(ax, 'Time (s)', 'FontSize', FontSize);
            else; set(ax, 'XTickLabel', []); end

            sig_str = '';
            if cdata.nEvents >= min_saccades && cdata.p_value < alpha_thresh
                sig_str = ' *';
            end
            title(ax, sprintf('%s  (n=%d)%s', col_labels{iC}, cdata.nEvents, sig_str), ...
                'FontSize', TitleFontSize, 'FontWeight', 'bold');
            box(ax, 'off');
        end
    end

    ax_label = axes('Position', [0 0 1 1], 'Visible', 'off');
    for iR = 1:2
        y_center = 1 - top_m - (iR-0.5)*row_h - (iR-1)*row_gap;
        text(ax_label, 0.02, y_center, row_labels{iR}, ...
            'FontSize', FontSize, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Rotation', 90, 'Units', 'normalized');
    end

    if strcmp(nd.pref_direction, 'CCW')
        bm_all_val = nd.all_CCW.burst_mod;
        bm_stat_val = nd.stat_CCW.burst_mod;
    elseif strcmp(nd.pref_direction, 'CW')
        bm_all_val = nd.all_CW.burst_mod;
        bm_stat_val = nd.stat_CW.burst_mod;
    else
        bm_all_val = NaN; bm_stat_val = NaN;
    end
    sgtitle(sprintf('%s — Pref: %s | BurstMod: %.2f (all) / %.2f (stat)', ...
        strrep(nd.neuron_name, '_', '\_'), nd.pref_direction, bm_all_val, bm_stat_val), ...
        'FontSize', TitleFontSize + 2, 'FontWeight', 'bold');
end

%% ========================================================================
%  LOCAL HELPER: population PETH figure (shared by workspace + full modes)
%  ========================================================================
function fh = make_population_figure(pop_pref_all, pop_pref_stat, ...
        pop_nonpref_all, pop_nonpref_stat, peth_centers, n_neurons, ...
        n_sessions, min_saccades, normalize_method, peth_window, ...
        FontSize, TickFontSize, TitleFontSize, doShow)

    % Compute mean +/- SEM
    pop_mean_arr = {mean(pop_pref_all, 1, 'omitnan'), mean(pop_nonpref_all, 1, 'omitnan'); ...
                    mean(pop_pref_stat, 1, 'omitnan'), mean(pop_nonpref_stat, 1, 'omitnan')};
    pop_sem_arr  = {std(pop_pref_all, 0, 1, 'omitnan') / sqrt(n_neurons), ...
                    std(pop_nonpref_all, 0, 1, 'omitnan') / sqrt(n_neurons); ...
                    std(pop_pref_stat, 0, 1, 'omitnan') / sqrt(n_neurons), ...
                    std(pop_nonpref_stat, 0, 1, 'omitnan') / sqrt(n_neurons)};

    if doShow
        fh = figure('Position', [50 50 1100 800]);
    else
        fh = figure('Position', [50 50 1100 800], 'Visible', 'off');
    end
    fh.UserData.created  = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fh.UserData.function = 'analyze_stationary_saccade_bursts_ver_1_1';
    fh.UserData.type     = 'population_PETH';

    pop_row_labels = {'All Saccades', 'Stationary Only'};
    pop_col_labels = {'Preferred', 'Non-Preferred'};

    % Distinct colors: row 1 (all) = saturated, row 2 (stationary) = deeper/darker
    line_colors    = {[0.20 0.50 0.80], [0.90 0.40 0.15]; ...     % all: bright blue, bright orange
                      [0.10 0.30 0.55], [0.65 0.25 0.05]};        % stat: dark blue, dark orange
    shade_clrs     = {[0.70 0.82 0.95], [0.97 0.80 0.65]; ...     % all: light blue, light orange
                      [0.35 0.50 0.70], [0.80 0.55 0.40]};        % stat: medium blue, medium orange

    % Shared y-axis
    all_ylim = [];
    for iR = 1:2
        for iC = 1:2
            all_ylim = [all_ylim, pop_mean_arr{iR,iC} + pop_sem_arr{iR,iC}, ...
                        pop_mean_arr{iR,iC} - pop_sem_arr{iR,iC}]; %#ok<AGROW>
        end
    end
    y_lo = min(all_ylim) - 0.5;
    y_hi = max(all_ylim) + 0.5;

    lm = 0.10; rm = 0.05; tm = 0.12; bm = 0.10;
    rg = 0.10; cg = 0.08;
    th = 1 - tm - bm;
    rh = (th - rg) / 2;
    tw = 1 - lm - rm;
    cw = (tw - cg) / 2;

    for iR = 1:2
        for iC = 1:2
            xp = lm + (iC-1)*(cw + cg);
            yp = 1 - tm - iR*rh - (iR-1)*rg;
            ax = axes('Position', [xp, yp, cw, rh]); %#ok<LAXES>
            hold(ax, 'on');

            m_tr = pop_mean_arr{iR, iC};
            s_tr = pop_sem_arr{iR, iC};

            fill_x = [peth_centers, fliplr(peth_centers)];
            fill_y = [m_tr + s_tr, fliplr(m_tr - s_tr)];
            fill(ax, fill_x, fill_y, shade_clrs{iR, iC}, ...
                'EdgeColor', 'none', 'FaceAlpha', 0.4);
            plot(ax, peth_centers, m_tr, '-', 'Color', line_colors{iR, iC}, ...
                'LineWidth', 4);

            xline(ax, 0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
            xlim(ax, peth_window);
            ylim(ax, [y_lo, y_hi]);
            ax.FontSize = TickFontSize;

            if iR == 2; xlabel(ax, 'Time (s)', 'FontSize', FontSize); end
            if iC == 1
                if strcmp(normalize_method, 'zscore')
                    ylabel(ax, 'FR (z-score)', 'FontSize', FontSize);
                else
                    ylabel(ax, 'FR (fold change)', 'FontSize', FontSize);
                end
            end

            % Panel title — no n here
            title(ax, sprintf('%s — %s', pop_row_labels{iR}, pop_col_labels{iC}), ...
                'FontSize', TitleFontSize, 'FontWeight', 'bold');
            box(ax, 'off');
        end
    end

    % Supertitle with n_neurons, n_sessions, and min saccades threshold
    sgtitle(sprintf('Population Saccade PETH — %s (%d neurons, %d sessions, >= %d saccades/dir)', ...
        normalize_method, n_neurons, n_sessions, min_saccades), ...
        'FontSize', TitleFontSize + 2, 'FontWeight', 'bold');
end

%% Version History
% 2026-04-04  v1.1  Refactored from v1.0:
%                    - min_saccades default changed from 5 to 8
%                    - "enough stationary" now requires BOTH CCW and CW
%                      directions to meet min_saccades (was OR)
%                    - Added pref_direction_list input (containers.Map)
%                      for user-supplied preferred directions
%                    - Added population PETH summary figure: 2x2 layout
%                      {All, Stationary} x {Preferred, Non-preferred},
%                      z-scored (or fold-change), mean +/- SEM shading
%                    - Stored preferred/non-preferred FR traces per neuron
%                      for population averaging
%                    - Added enough_stat_both column to summary table
%                    - Population PETH data stored in params struct
% 2026-04-03  v1.0  Initial version. Per-trial Wilcoxon signed-rank test.
%                    Burst modulation index. Per-neuron diagnostic figures.
%                    NOTE: Using Jeff's saccade times from
%                    saccades-edited.mat. Chenhao's revised onset times
%                    pending.
