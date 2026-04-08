function [batch_summary, all_results, params] = analyze_stationary_saccade_bursts_ver_1_0(varargin)
% analyze_stationary_saccade_bursts_ver_1_0  Population analysis of
% saccade-evoked bursting during platform-stationary vs. moving periods.
%
%   For each neuron in the input list (default: 104 NPH neurons), this
%   function:
%     1. Loads spikes, saccade times, and stationary period annotations
%     2. Separates saccades into moving vs. stationary
%     3. Computes saccade-aligned FR PETHs for each condition x direction
%     4. Determines the preferred burst direction (from all-saccade data)
%     5. Quantifies burst modulation and tests significance via Wilcoxon
%        signed-rank test (per-trial burst window FR vs. baseline FR)
%     6. Generates a diagnostic figure per neuron (optional)
%     7. Saves all results + figures to a single output folder
%
%   INPUTS:
%       varargin    - name-value pairs (see defaults below)
%
%   KEY OPTIONS:
%       cell_list           - cell array of full t-file paths (default:
%                             loads NPH_all_tfilelist from tfilelist_path)
%       tfilelist_path      - path to NPH all tfilelist.mat
%       data_root           - root data directory (used if cell_list is
%                             empty and tfilelist_path is empty: scans all
%                             session folders)
%       save_root           - analysis output path
%       min_saccades        - minimum stationary saccades per direction to
%                             include in analysis (default: 5)
%       burst_window        - [start end] in s relative to saccade (default: [-0.05 0.05])
%       baseline_window     - [start end] in s relative to saccade (default: [-0.2 -0.1])
%       alpha_threshold     - significance level (default: 0.05)
%       peth_window         - full PETH window (default: [-0.2 0.2])
%       peth_dt             - PETH bin width in s (default: 0.005)
%       peth_smooth_win     - smoothdata gaussian window in bins (default: 5)
%       doMakeFigs          - generate per-neuron figures (default: 1)
%       doSavePNG, doSaveFIG
%
%   OUTPUTS:
%       batch_summary   - MATLAB table with one row per neuron:
%                         neuron_name, SSN, nCCW_all, nCW_all,
%                         nCCW_stat, nCW_stat, nCCW_move, nCW_move,
%                         pref_direction, burst_mod_all, burst_mod_move,
%                         burst_mod_stat, p_all, p_move, p_stat,
%                         sig_all, sig_move, sig_stat,
%                         stationary_dur_s, n_stat_periods, status
%       all_results     - cell array of per-neuron result structs
%       params          - struct of parameter values used
%
%   USAGE:
%       % Run on all 104 NPH neurons:
%       [summary, results, par] = analyze_stationary_saccade_bursts_ver_1_0();
%
%       % Run on a custom list:
%       my_list = NPH_all_tfilelist(1:10);
%       [summary, results, par] = analyze_stationary_saccade_bursts_ver_1_0('cell_list', my_list);
%
%   DEPENDENCIES:
%       LoadSpikesQuick, SelectTS, restrict, pushdir, popdir,
%       HD_GetSSN, process_varargin, tsd, nearest_idx3, MakeQfromS

%% ========================================================================
%  DEFAULTS
%  ========================================================================
cell_list           = {};       % if empty, loads from tfilelist_path
tfilelist_path      = 'C:\Jeff\U01\datatouse\NPH all tfilelist.mat';
data_root           = 'C:\Jeff\U01\datatouse';
save_root           = 'C:\Jeff\U01\datatouse\___analysis_output';

% Analysis windows
burst_window        = [-0.05 0.05];     % s, peri-saccade burst measurement
baseline_window     = [-0.2 -0.1];      % s, pre-saccade baseline
alpha_threshold     = 0.05;             % significance level

% Saccade PETH parameters
peth_window         = [-0.2 0.2];
peth_dt             = 0.005;
peth_smooth_win     = 5;
peth_excessBounds   = 1;

% Stationary saccade filter
min_saccades        = 5;

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
doMakeFigs  = 1;
doSavePNG   = 1;
doSaveFIG   = 1;
doShow      = 0;
verbose     = 1;

%% Process varargin
process_varargin(varargin);

%% ========================================================================
%  PARAMS (analysis-relevant free parameters only)
%  ========================================================================
params = struct();
params.function_name    = 'analyze_stationary_saccade_bursts';
params.version          = '1.0';
params.date_run         = datestr(now, 'yyyy-mm-dd HH:MM:SS');
params.burst_window     = burst_window;
params.baseline_window  = baseline_window;
params.alpha_threshold  = alpha_threshold;
params.peth_window      = peth_window;
params.peth_dt          = peth_dt;
params.peth_smooth_win  = peth_smooth_win;
params.peth_excessBounds = peth_excessBounds;
params.min_saccades     = min_saccades;

%% ========================================================================
%  BUILD CELL LIST
%  ========================================================================
if isempty(cell_list)
    if exist(tfilelist_path, 'file')
        tfl = load(tfilelist_path, 'NPH_all_tfilelist');
        cell_list = tfl.NPH_all_tfilelist(:);  % ensure column
        if verbose; fprintf('Loaded %d neurons from tfilelist.\n', length(cell_list)); end
    else
        % Scan data_root for all .t files
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
sum_status          = cell(nNeurons, 1);

all_results = cell(nNeurons, 1);

%% ========================================================================
%  MAIN LOOP
%  ========================================================================
peth_bins    = linspace(peth_window(1), peth_window(2), diff(peth_window)/peth_dt + 1);
peth_centers = peth_bins(1:end-1) + diff(peth_bins) / 2;

for iN = 1:nNeurons
    tfile_path = cell_list{iN};
    [sess_path, neuron_name] = fileparts(tfile_path);

    if verbose
        fprintf('\n=== %d/%d: %s ===\n', iN, nNeurons, neuron_name);
    end

    % Initialize result struct
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
        [S_all, spike_cell_names, spk_params] = LoadSpikesQuick();

        % Match cell
        cell_idx = find_cell_index(spike_cell_names, neuron_name, SSN);
        if isempty(cell_idx)
            popdir;
            sum_status{iN} = 'ERROR: cell not found in spike files';
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
            if verbose; fprintf('  SKIP: no saccade file\n'); end
            all_results{iN} = nd;
            continue
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
            if verbose; fprintf('  SKIP: no stationary times file\n'); end
            all_results{iN} = nd;
            continue
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

        nd.nCCW_all  = nT_all;     nd.nCW_all  = nN_all;
        nd.nCCW_stat = nT_stat;    nd.nCW_stat = nN_stat;
        nd.nCCW_move = nT_move;    nd.nCW_move = nN_move;
        nd.stat_dur  = stat_dur_total;
        nd.n_stat_periods = length(STtstart);

        sum_nCCW_all(iN)  = nT_all;    sum_nCW_all(iN)  = nN_all;
        sum_nCCW_stat(iN) = nT_stat;   sum_nCW_stat(iN) = nN_stat;
        sum_nCCW_move(iN) = nT_move;   sum_nCW_move(iN) = nN_move;
        sum_stat_dur(iN)  = stat_dur_total;
        sum_n_stat_periods(iN) = length(STtstart);

        if verbose
            fprintf('  CCW: %d all | %d move | %d stat\n', nT_all, nT_move, nT_stat);
            fprintf('  CW:  %d all | %d move | %d stat\n', nN_all, nN_move, nN_stat);
            fprintf('  Stationary: %d periods, %.1f s total\n', length(STtstart), stat_dur_total);
        end

        %% ---- Compute PETHs and burst metrics for all 6 conditions ----
        % Conditions: {all, move, stat} x {CCW, CW}
        conditions = struct( ...
            'label',     {'all_CCW', 'all_CW', 'move_CCW', 'move_CW', 'stat_CCW', 'stat_CW'}, ...
            'sacc_times', {temporal_all, nasal_all, temporal_move, nasal_move, temporal_stat, nasal_stat}, ...
            'nEvents',   {nT_all, nN_all, nT_move, nN_move, nT_stat, nN_stat} ...
        );

        for iCond = 1:6
            cond = conditions(iCond);

            if cond.nEvents < 1
                nd.(cond.label).fr_smooth   = nan(size(peth_centers));
                nd.(cond.label).burst_mod   = NaN;
                nd.(cond.label).p_value     = NaN;
                nd.(cond.label).is_sig      = false;
                nd.(cond.label).mean_burst_fr    = NaN;
                nd.(cond.label).mean_baseline_fr = NaN;
                nd.(cond.label).nEvents     = 0;
                nd.(cond.label).spk_rel     = [];
                nd.(cond.label).evt_idx     = [];
                continue
            end

            % Raster data
            [spk_rel, evt_idx] = compute_raster_data(myCell, cond.sacc_times, ...
                peth_window, peth_excessBounds);

            % Smoothed FR PETH
            [m_hist, ~] = histcounts(spk_rel, peth_bins);
            fr_peth = m_hist / peth_dt / max(cond.nEvents, 1);
            fr_smooth = smoothdata(movmean(fr_peth, peth_smooth_win), 'gaussian', peth_smooth_win);

            % Per-trial burst and baseline FR for significance test
            [trial_burst_fr, trial_baseline_fr] = compute_per_trial_fr( ...
                myCell, cond.sacc_times, burst_window, baseline_window);

            % Burst modulation index
            mean_burst_fr    = mean(trial_burst_fr, 'omitnan');
            mean_baseline_fr = mean(trial_baseline_fr, 'omitnan');
            if mean_baseline_fr > 0
                burst_mod_val = (mean_burst_fr - mean_baseline_fr) / mean_baseline_fr;
            else
                burst_mod_val = NaN;
            end

            % Wilcoxon signed-rank test (paired: burst vs baseline per trial)
            valid_trials = ~isnan(trial_burst_fr) & ~isnan(trial_baseline_fr);
            if sum(valid_trials) >= 3
                [p_val, ~] = signrank(trial_burst_fr(valid_trials), ...
                    trial_baseline_fr(valid_trials));
            else
                p_val = NaN;
            end

            % Store
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

        %% ---- Determine preferred direction (from ALL saccades) ----
        mod_CCW = nd.all_CCW.burst_mod;
        mod_CW  = nd.all_CW.burst_mod;

        if isnan(mod_CCW) && isnan(mod_CW)
            pref_dir = 'none';
        elseif isnan(mod_CW) || mod_CCW >= mod_CW
            pref_dir = 'CCW';
        else
            pref_dir = 'CW';
        end
        nd.pref_direction = pref_dir;
        sum_pref_dir{iN} = pref_dir;

        % Extract preferred-direction metrics for summary
        if strcmp(pref_dir, 'CCW')
            sum_burst_mod_all(iN)  = nd.all_CCW.burst_mod;
            sum_burst_mod_move(iN) = nd.move_CCW.burst_mod;
            sum_burst_mod_stat(iN) = nd.stat_CCW.burst_mod;
            sum_p_all(iN)          = nd.all_CCW.p_value;
            sum_p_move(iN)         = nd.move_CCW.p_value;
            sum_p_stat(iN)         = nd.stat_CCW.p_value;
        elseif strcmp(pref_dir, 'CW')
            sum_burst_mod_all(iN)  = nd.all_CW.burst_mod;
            sum_burst_mod_move(iN) = nd.move_CW.burst_mod;
            sum_burst_mod_stat(iN) = nd.stat_CW.burst_mod;
            sum_p_all(iN)          = nd.all_CW.p_value;
            sum_p_move(iN)         = nd.move_CW.p_value;
            sum_p_stat(iN)         = nd.stat_CW.p_value;
        end

        sum_sig_all(iN)  = sum_p_all(iN) < alpha_threshold;
        sum_sig_move(iN) = sum_p_move(iN) < alpha_threshold;
        sum_sig_stat(iN) = sum_p_stat(iN) < alpha_threshold;

        % Check if enough stationary saccades
        enough_stat = (nT_stat >= min_saccades) || (nN_stat >= min_saccades);
        nd.enough_stationary = enough_stat;

        nd.peth_centers = peth_centers;

        if verbose
            fprintf('  Preferred: %s | Burst mod (all/move/stat): %.2f / %.2f / %.2f\n', ...
                pref_dir, sum_burst_mod_all(iN), sum_burst_mod_move(iN), sum_burst_mod_stat(iN));
            fprintf('  p-values (all/move/stat): %.4f / %.4f / %.4f\n', ...
                sum_p_all(iN), sum_p_move(iN), sum_p_stat(iN));
        end

        %% ---- Generate figure ----
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
    sum_stat_dur, sum_n_stat_periods, sum_status, ...
    'VariableNames', { ...
    'neuron_name', 'SSN', ...
    'nCCW_all', 'nCW_all', 'nCCW_stat', 'nCW_stat', ...
    'nCCW_move', 'nCW_move', ...
    'pref_direction', 'burst_mod_all', 'burst_mod_move', 'burst_mod_stat', ...
    'p_all', 'p_move', 'p_stat', ...
    'sig_all', 'sig_move', 'sig_stat', ...
    'stationary_dur_s', 'n_stat_periods', 'status'});

%% ========================================================================
%  PRINT SUMMARY
%  ========================================================================
n_ok = sum(strcmp(sum_status, 'OK'));
n_enough_stat = sum(sum_nCCW_stat >= min_saccades | sum_nCW_stat >= min_saccades);
n_sig_all  = sum(sum_sig_all);
n_sig_move = sum(sum_sig_move);
n_sig_stat = sum(sum_sig_stat & (sum_nCCW_stat >= min_saccades | sum_nCW_stat >= min_saccades));

fprintf('\n\n====== POPULATION SUMMARY ======\n');
fprintf('  Total neurons:          %d\n', nNeurons);
fprintf('  Successfully processed: %d\n', n_ok);
fprintf('  With >= %d stat saccades: %d\n', min_saccades, n_enough_stat);
fprintf('  Sig burst (all sacc):   %d / %d (%.0f%%)\n', n_sig_all, n_ok, 100*n_sig_all/max(n_ok,1));
fprintf('  Sig burst (moving):     %d / %d (%.0f%%)\n', n_sig_move, n_ok, 100*n_sig_move/max(n_ok,1));
fprintf('  Sig burst (stationary): %d / %d with enough saccades (%.0f%%)\n', ...
    n_sig_stat, n_enough_stat, 100*n_sig_stat/max(n_enough_stat,1));
fprintf('================================\n');

%% ========================================================================
%  SAVE
%  ========================================================================
save(fullfile(out_folder, 'stationary_burst_analysis.mat'), ...
    'batch_summary', 'all_results', 'params', '-v7.3');
if verbose; fprintf('\nResults saved to: %s\n', out_folder); end

% Also save summary as CSV for quick inspection
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
% For each saccade event, count spikes in the burst window and baseline
% window, convert to firing rate. Returns vectors of length nEvents.

    nEvents = length(saccade_times);
    trial_burst_fr    = nan(nEvents, 1);
    trial_baseline_fr = nan(nEvents, 1);

    burst_dur    = diff(burst_win);
    baseline_dur = diff(baseline_win);

    for iE = 1:nEvents
        t_sacc = saccade_times(iE);

        % Burst window spikes
        S_burst = restrict(myCell, t_sacc + burst_win(1), t_sacc + burst_win(2));
        trial_burst_fr(iE) = length(S_burst.t{1}) / burst_dur;

        % Baseline window spikes
        S_base = restrict(myCell, t_sacc + baseline_win(1), t_sacc + baseline_win(2));
        trial_baseline_fr(iE) = length(S_base.t{1}) / baseline_dur;
    end
end

%% ========================================================================
%  LOCAL HELPER: find cell index in spike_cell_names
%  ========================================================================
function cell_idx = find_cell_index(spike_cell_names, neuron_name, SSN)
    cell_idx = [];
    for cc = 1:length(spike_cell_names)
        if strcmp(spike_cell_names{cc}, neuron_name)
            cell_idx = cc;
            return
        end
    end
    % Try matching without session prefix
    for cc = 1:length(spike_cell_names)
        if contains(spike_cell_names{cc}, strrep(neuron_name, [SSN '-'], ''))
            cell_idx = cc;
            return
        end
    end
end

%% ========================================================================
%  LOCAL HELPER: generate diagnostic figure
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

    % Store metadata
    fh.UserData.created   = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fh.UserData.function  = 'analyze_stationary_saccade_bursts_ver_1_0';
    fh.UserData.neuron    = nd.neuron_name;

    % Layout: 2 rows x 2 cols
    row_labels = {'All Saccades', 'Stationary Only'};
    col_labels = {'CCW (Temporal)', 'CW (Nasal)'};

    cond_keys = {'all_CCW', 'all_CW'; 'stat_CCW', 'stat_CW'};
    fr_colors = {color_all_CCW, color_all_CW; color_stat_CCW, color_stat_CW};

    % Compute shared FR y-axis
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

    % Axes positioning
    left_margin   = 0.10;  right_margin  = 0.08;
    top_margin    = 0.12;  bottom_margin = 0.10;
    row_gap = 0.08;        col_gap = 0.06;

    total_height = 1 - top_margin - bottom_margin;
    row_height   = (total_height - row_gap) / 2;
    total_width  = 1 - left_margin - right_margin;
    col_width    = (total_width - col_gap) / 2;

    for iR = 1:2
        for iC = 1:2
            x_pos = left_margin + (iC - 1) * (col_width + col_gap);
            y_pos = 1 - top_margin - iR * row_height - (iR - 1) * row_gap;

            ax = axes('Position', [x_pos, y_pos, col_width, row_height]); %#ok<LAXES>
            hold(ax, 'on');

            key = cond_keys{iR, iC};
            cdata = nd.(key);

            if cdata.nEvents >= min_saccades && ~isempty(cdata.spk_rel)
                % Raster
                yyaxis(ax, 'left');
                plot(ax, cdata.spk_rel, cdata.evt_idx + 0.5, '.', ...
                    'Color', rasterColor, 'MarkerSize', MarkerSize_raster);
                ylim(ax, [0 cdata.nEvents + 1]);
                ax.YAxis(1).Color = [0 0 0];
                if iC == 1; ylabel(ax, 'Trial', 'FontSize', FontSize);
                else; set(ax.YAxis(1), 'TickLabels', []); end

                % FR overlay
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
                text(ax, 0.5, 0.5, 'n = 0', ...
                    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                    'FontSize', FontSize, 'Color', [0.6 0.6 0.6]);
            end

            xline(ax, 0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
            xlim(ax, peth_window);
            ax.XAxis.FontSize = TickFontSize;

            if iR == 2; xlabel(ax, 'Time (s)', 'FontSize', FontSize);
            else; set(ax, 'XTickLabel', []); end

            % Title with significance marker
            sig_str = '';
            if cdata.nEvents >= min_saccades && cdata.p_value < alpha_thresh
                sig_str = ' *';
            end
            title_str = sprintf('%s  (n=%d)%s', col_labels{iC}, cdata.nEvents, sig_str);
            title(ax, title_str, 'FontSize', TitleFontSize, 'FontWeight', 'bold');

            box(ax, 'off');
        end
    end

    % Row labels
    ax_label = axes('Position', [0 0 1 1], 'Visible', 'off');
    for iR = 1:2
        y_center = 1 - top_margin - (iR - 0.5) * row_height - (iR - 1) * row_gap;
        text(ax_label, 0.02, y_center, row_labels{iR}, ...
            'FontSize', FontSize, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Rotation', 90, 'Units', 'normalized');
    end

    % Supertitle
    if strcmp(nd.pref_direction, 'CCW')
        bm_all_val = nd.all_CCW.burst_mod;
        bm_stat_val = nd.stat_CCW.burst_mod;
    elseif strcmp(nd.pref_direction, 'CW')
        bm_all_val = nd.all_CW.burst_mod;
        bm_stat_val = nd.stat_CW.burst_mod;
    else
        bm_all_val = NaN; bm_stat_val = NaN;
    end
    sgtitle_str = sprintf('%s — Pref: %s | BurstMod: %.2f (all) / %.2f (stat)', ...
        strrep(nd.neuron_name, '_', '\_'), nd.pref_direction, bm_all_val, bm_stat_val);
    sgtitle(sgtitle_str, 'FontSize', TitleFontSize + 2, 'FontWeight', 'bold');
end

%% Version History
% 2026-04-03  v1.0  Initial version. Population analysis of saccade-evoked
%                    bursting during stationary vs. moving periods across
%                    all 104 NPH neurons. Per-trial Wilcoxon signed-rank
%                    test for significance (burst window vs. baseline).
%                    Burst modulation index = (burst-baseline)/baseline.
%                    Preferred direction determined from all-saccade
%                    condition. Outputs: summary table, per-neuron results
%                    struct, figures. NOTE: Using Jeff's saccade times from
%                    saccades-edited.mat. Chenhao's revised onset times
%                    pending.
