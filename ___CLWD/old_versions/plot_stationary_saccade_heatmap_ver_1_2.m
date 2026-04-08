function [Z, params] = plot_stationary_saccade_heatmap_ver_1_2(varargin)
% plot_stationary_saccade_heatmap_ver_1_2  v1.2
% 4-panel population heatmap of saccade PETHs: VOR vs Spontaneous,
% Preferred vs Non-Preferred.
%
%   Top row:    VOR Preferred     | VOR Non-Preferred
%   Bottom row: Spont Preferred   | Spont Non-Preferred
%
%   All panels share the same row ordering (sorted by VOR preferred
%   peak bin), so visual comparison across conditions is straightforward.
%
%   Uses Chenhao's saccade onset times (saccade_ONSET_all_unit.mat).
%   Preferred direction determined from VOR saccade burst modulation.
%   Neurons must have >= min_saccades spontaneous saccades in BOTH
%   directions to qualify.
%
%   INPUTS (name-value pairs via varargin):
%     tfilelist_path    - path to NPH_all_tfilelist.mat
%     min_saccades      - min stationary saccades per dir (default: 10)
%     precomputed_Z     - Z struct from prior run (workspace mode)
%     peth_window       - [start end] in s (default: [-0.2 0.2])
%     peth_dt           - bin width in s (default: 0.005)
%     peth_smooth_win   - smoothing window in bins (default: 5)
%     burst_window      - for pref dir detection (default: [-0.05 0.05])
%     baseline_window   - for pref dir detection (default: [-0.2 -0.1])
%     doSavePNG/FIG/EPS - save flags (default: 1)
%     doShow            - display figure (default: 1)
%     save_dir          - output root (default: ___analysis_output)
%
%   OUTPUTS:
%     Z      - struct with FR matrices, sort indices, bin centers
%     params - struct of all free parameter values
%
%   DEPENDENCIES:
%     LoadSpikesQuick, SelectTS, restrict, pushdir, popdir, HD_GetSSN,
%     process_varargin, load_chenhao_saccades, save_fig_trio, find_cell_index

%% =====================================================================
%  DEFAULTS
%  =====================================================================
tfilelist_path    = 'C:\Jeff\U01\datatouse\NPH all tfilelist.mat';
save_dir          = 'C:\Jeff\U01\datatouse\___analysis_output';
precomputed_Z     = [];

% PETH parameters
peth_window       = [-0.2 0.2];
peth_dt           = 0.005;
peth_smooth_win   = 5;
peth_excessBounds = 1;

% Burst modulation windows (for preferred direction detection)
burst_window      = [-0.05 0.05];
baseline_window   = [-0.2 -0.1];

% Thresholds
min_saccades      = 10;

% Figure
FontSize          = 28;
TitleFontSize     = 24;
TickFontSize      = 18;
LineWidth_zero    = 2;

% Output control
doSavePNG         = 1;
doSaveFIG         = 1;
doSaveEPS         = 1;
doShow            = 1;
verbose           = 1;

%% Process varargin
process_varargin(varargin);

% Manual extraction for struct params (process_varargin can silently fail on structs)
for iV = 1:2:length(varargin)
    if ischar(varargin{iV}) && strcmp(varargin{iV}, 'precomputed_Z')
        precomputed_Z = varargin{iV+1};
    end
end

%% =====================================================================
%  PARAMS
%  =====================================================================
params = struct();
params.function_name    = 'plot_stationary_saccade_heatmap';
params.version          = '1.2';
params.date_run         = datestr(now, 'yyyy-mm-dd HH:MM:SS');
params.peth_window      = peth_window;
params.peth_dt          = peth_dt;
params.peth_smooth_win  = peth_smooth_win;
params.burst_window     = burst_window;
params.baseline_window  = baseline_window;
params.min_saccades     = min_saccades;
params.saccade_source   = 'Chenhao onset times (saccade_ONSET_all_unit.mat)';

%% =====================================================================
%  SETUP
%  =====================================================================
peth_bins    = linspace(peth_window(1), peth_window(2), diff(peth_window)/peth_dt + 1);
peth_centers = peth_bins(1:end-1) + diff(peth_bins) / 2;
nBins        = length(peth_centers);
burst_mask   = peth_centers >= burst_window(1) & peth_centers <= burst_window(2);
base_mask    = peth_centers >= baseline_window(1) & peth_centers <= baseline_window(2);

%% =====================================================================
%  WORKSPACE MODE — skip recomputation if Z struct provided
%  =====================================================================
if ~isempty(precomputed_Z) && isfield(precomputed_Z, 'FR_vor_pref_smooth')
    fprintf('=== WORKSPACE MODE (Z struct with %d neurons) ===\n', precomputed_Z.n_neurons);
    Z = precomputed_Z;

    % Allow re-thresholding from workspace
    if isfield(Z, 'all_FR_stat_pref') && min_saccades ~= Z.min_saccades
        fprintf('  Re-thresholding: %d -> %d min saccades\n', Z.min_saccades, min_saccades);
        new_keep = Z.all_n_stat_CCW >= min_saccades & Z.all_n_stat_CW >= min_saccades & Z.all_qualifies;

        Z.FR_stat_pref_smooth    = smoothdata(Z.all_FR_stat_pref(new_keep, :), 2, 'gaussian', peth_smooth_win);
        Z.FR_stat_nonpref_smooth = smoothdata(Z.all_FR_stat_nonpref(new_keep, :), 2, 'gaussian', peth_smooth_win);
        Z.FR_vor_pref_smooth     = smoothdata(Z.all_FR_vor_pref(new_keep, :), 2, 'gaussian', peth_smooth_win);
        Z.FR_vor_nonpref_smooth  = smoothdata(Z.all_FR_vor_nonpref(new_keep, :), 2, 'gaussian', peth_smooth_win);
        Z.neuron_names           = Z.all_neuron_names(new_keep);
        Z.n_neurons              = sum(new_keep);
        Z.min_saccades           = min_saccades;

        [~, peak_bin_re] = max(Z.FR_vor_pref_smooth, [], 2);
        [~, Z.sort_by_vor_pref] = sort(peak_bin_re);
    end

    % Output folder
    out_folder = make_output_folder(save_dir, doSavePNG, doSaveFIG, doSaveEPS);
    params.output_folder = out_folder;
    peth_centers = Z.binCenters;

    % Plot
    plot_4panel_heatmap(Z.FR_vor_pref_smooth, Z.FR_vor_nonpref_smooth, ...
        Z.FR_stat_pref_smooth, Z.FR_stat_nonpref_smooth, ...
        Z.sort_by_vor_pref, peth_centers, Z.n_neurons, Z.n_sessions, ...
        min_saccades, FontSize, TitleFontSize, TickFontSize, LineWidth_zero, ...
        doShow, doSavePNG, doSaveFIG, doSaveEPS, out_folder);

    if doSavePNG || doSaveFIG || doSaveEPS
        save(fullfile(out_folder, 'heatmap_data.mat'), 'Z', 'params', '-v7.3');
        fprintf('Saved to: %s\n', out_folder);
    end
    return
end

%% =====================================================================
%  FULL RUN MODE
%  =====================================================================

% Load neuron list
tfl = load(tfilelist_path, 'NPH_all_tfilelist');
cell_list = tfl.NPH_all_tfilelist(:);
nNeurons = length(cell_list);
if verbose; fprintf('Loaded %d neurons from tfilelist.\n', nNeurons); end

% Output folder
out_folder = make_output_folder(save_dir, doSavePNG, doSaveFIG, doSaveEPS);
params.output_folder = out_folder;

%% =====================================================================
%  PRE-ALLOCATE
%  =====================================================================
FR_stat_pref     = nan(nNeurons, nBins);
FR_stat_nonpref  = nan(nNeurons, nBins);
FR_vor_pref      = nan(nNeurons, nBins);
FR_vor_nonpref   = nan(nNeurons, nBins);
neuron_names     = cell(nNeurons, 1);
session_names    = cell(nNeurons, 1);
pref_dirs        = cell(nNeurons, 1);
n_stat_CCW       = nan(nNeurons, 1);
n_stat_CW        = nan(nNeurons, 1);
n_vor_CCW        = nan(nNeurons, 1);
n_vor_CW         = nan(nNeurons, 1);
qualifies        = false(nNeurons, 1);
prev_SSN         = '';

%% =====================================================================
%  MAIN LOOP
%  =====================================================================
for iN = 1:nNeurons
    tfile_path = cell_list{iN};
    [sess_path, tfile_name] = fileparts(tfile_path);
    neuron_names{iN} = tfile_name;

    if verbose && mod(iN, 10) == 0
        fprintf('  %d/%d ...\n', iN, nNeurons);
    end

    try
        pushdir(sess_path);
        SSN = HD_GetSSN;
        session_names{iN} = SSN;

        %% Load session data (only reload if new session)
        if ~strcmp(SSN, prev_SSN)
            [S_all, spike_cell_names] = LoadSpikesQuick();

            % Load Chenhao saccade times
            [temporal_all, nasal_all] = load_chenhao_saccades(SSN);
            if isempty(temporal_all) && isempty(nasal_all)
                popdir; prev_SSN = SSN; continue
            end

            % Load stationary times
            stat_file = strcat(SSN, '-AHV_StationaryTimes.mat');
            if ~exist(stat_file, 'file'); popdir; prev_SSN = SSN; continue; end
            stat_data = load(stat_file, 'STtstart', 'STtend');
            STtstart = stat_data.STtstart(:);
            STtend   = stat_data.STtend(:);

            % Split VOR vs stationary saccades
            t_stat_mask = false(size(temporal_all));
            n_stat_mask = false(size(nasal_all));
            for iST = 1:length(STtstart)
                t_stat_mask = t_stat_mask | (temporal_all > STtstart(iST) & temporal_all < STtend(iST));
                n_stat_mask = n_stat_mask | (nasal_all > STtstart(iST) & nasal_all < STtend(iST));
            end
            sess_temporal_stat = temporal_all(t_stat_mask);
            sess_nasal_stat    = nasal_all(n_stat_mask);
            sess_temporal_move = temporal_all(~t_stat_mask);
            sess_nasal_move    = nasal_all(~n_stat_mask);

            prev_SSN = SSN;
        end

        % Find this neuron in the spike data
        cell_idx_val = find_cell_index(spike_cell_names, tfile_name, SSN);
        if isempty(cell_idx_val); popdir; continue; end
        myCell = SelectTS([], S_all, cell_idx_val);

        n_stat_CCW(iN) = length(sess_temporal_stat);
        n_stat_CW(iN)  = length(sess_nasal_stat);
        n_vor_CCW(iN)  = length(sess_temporal_move);
        n_vor_CW(iN)   = length(sess_nasal_move);

        %% Compute 4 PETHs
        fr_vor_CCW_trace  = compute_peth(myCell, sess_temporal_move, peth_bins, peth_dt, peth_smooth_win, peth_excessBounds);
        fr_vor_CW_trace   = compute_peth(myCell, sess_nasal_move, peth_bins, peth_dt, peth_smooth_win, peth_excessBounds);
        fr_stat_CCW_trace = compute_peth(myCell, sess_temporal_stat, peth_bins, peth_dt, peth_smooth_win, peth_excessBounds);
        fr_stat_CW_trace  = compute_peth(myCell, sess_nasal_stat, peth_bins, peth_dt, peth_smooth_win, peth_excessBounds);

        %% Determine preferred direction from VOR burst modulation
        bm_CCW = compute_burst_mod_value(fr_vor_CCW_trace, base_mask, burst_mask);
        bm_CW  = compute_burst_mod_value(fr_vor_CW_trace, base_mask, burst_mask);

        if isnan(bm_CCW) && isnan(bm_CW)
            pref_dir_str = 'none';
        elseif isnan(bm_CW) || bm_CCW >= bm_CW
            pref_dir_str = 'CCW';
        else
            pref_dir_str = 'CW';
        end
        pref_dirs{iN} = pref_dir_str;

        if strcmp(pref_dir_str, 'CCW')
            FR_stat_pref(iN,:)    = fr_stat_CCW_trace;
            FR_stat_nonpref(iN,:) = fr_stat_CW_trace;
            FR_vor_pref(iN,:)     = fr_vor_CCW_trace;
            FR_vor_nonpref(iN,:)  = fr_vor_CW_trace;
        else
            FR_stat_pref(iN,:)    = fr_stat_CW_trace;
            FR_stat_nonpref(iN,:) = fr_stat_CCW_trace;
            FR_vor_pref(iN,:)     = fr_vor_CW_trace;
            FR_vor_nonpref(iN,:)  = fr_vor_CCW_trace;
        end

        %% Check if qualifies
        qualifies(iN) = n_stat_CCW(iN) >= min_saccades && n_stat_CW(iN) >= min_saccades ...
            && ~strcmp(pref_dir_str, 'none');

        popdir;
    catch ME
        if verbose; fprintf('  ERROR %s: %s\n', tfile_name, ME.message); end
        try popdir; catch; end
    end
end

%% =====================================================================
%  FILTER TO QUALIFYING NEURONS
%  =====================================================================
keep_idx = find(qualifies);
nQual = length(keep_idx);
n_sessions = length(unique(session_names(keep_idx)));

fprintf('\n=== HEATMAP SUMMARY ===\n');
fprintf('  Total neurons:    %d\n', nNeurons);
fprintf('  Qualifying:       %d (>= %d stationary saccades both dirs)\n', nQual, min_saccades);
fprintf('  Sessions:         %d\n', n_sessions);
fprintf('=======================\n');

if nQual < 2
    warning('Not enough qualifying neurons (%d).', nQual);
    Z = struct(); return
end

% Extract qualifying data
qual_stat_pref    = FR_stat_pref(keep_idx, :);
qual_stat_nonpref = FR_stat_nonpref(keep_idx, :);
qual_vor_pref     = FR_vor_pref(keep_idx, :);
qual_vor_nonpref  = FR_vor_nonpref(keep_idx, :);

%% =====================================================================
%  SMOOTH + SORT
%  =====================================================================
qual_stat_pref_smooth    = smoothdata(qual_stat_pref, 2, 'gaussian', peth_smooth_win);
qual_stat_nonpref_smooth = smoothdata(qual_stat_nonpref, 2, 'gaussian', peth_smooth_win);
qual_vor_pref_smooth     = smoothdata(qual_vor_pref, 2, 'gaussian', peth_smooth_win);
qual_vor_nonpref_smooth  = smoothdata(qual_vor_nonpref, 2, 'gaussian', peth_smooth_win);

% Sort by peak bin of VOR preferred (consistent row order across all 4 panels)
[~, peak_bin_vor_pref] = max(qual_vor_pref_smooth, [], 2);
[~, sort_by_vor_pref]  = sort(peak_bin_vor_pref);

%% =====================================================================
%  BUILD Z STRUCT
%  =====================================================================
Z = struct();
Z.FR_stat_pref           = qual_stat_pref;
Z.FR_stat_nonpref        = qual_stat_nonpref;
Z.FR_stat_pref_smooth    = qual_stat_pref_smooth;
Z.FR_stat_nonpref_smooth = qual_stat_nonpref_smooth;
Z.FR_vor_pref            = qual_vor_pref;
Z.FR_vor_nonpref         = qual_vor_nonpref;
Z.FR_vor_pref_smooth     = qual_vor_pref_smooth;
Z.FR_vor_nonpref_smooth  = qual_vor_nonpref_smooth;
Z.sort_by_vor_pref       = sort_by_vor_pref;
Z.binCenters             = peth_centers;
Z.peth_bins              = peth_bins;
Z.neuron_names           = neuron_names(keep_idx);
Z.session_names          = session_names(keep_idx);
Z.pref_dirs              = pref_dirs(keep_idx);
Z.n_stat_CCW             = n_stat_CCW(keep_idx);
Z.n_stat_CW              = n_stat_CW(keep_idx);
Z.n_vor_CCW              = n_vor_CCW(keep_idx);
Z.n_vor_CW               = n_vor_CW(keep_idx);
Z.n_neurons              = nQual;
Z.n_sessions             = n_sessions;
Z.min_saccades           = min_saccades;

% All neuron data (unfiltered) for re-thresholding in workspace mode
Z.all_FR_stat_pref       = FR_stat_pref;
Z.all_FR_stat_nonpref    = FR_stat_nonpref;
Z.all_FR_vor_pref        = FR_vor_pref;
Z.all_FR_vor_nonpref     = FR_vor_nonpref;
Z.all_qualifies          = qualifies;
Z.all_neuron_names       = neuron_names;
Z.all_n_stat_CCW         = n_stat_CCW;
Z.all_n_stat_CW          = n_stat_CW;

%% =====================================================================
%  PLOT + SAVE
%  =====================================================================
plot_4panel_heatmap(qual_vor_pref_smooth, qual_vor_nonpref_smooth, ...
    qual_stat_pref_smooth, qual_stat_nonpref_smooth, ...
    sort_by_vor_pref, peth_centers, nQual, n_sessions, min_saccades, ...
    FontSize, TitleFontSize, TickFontSize, LineWidth_zero, ...
    doShow, doSavePNG, doSaveFIG, doSaveEPS, out_folder);

if doSavePNG || doSaveFIG || doSaveEPS
    save(fullfile(out_folder, 'heatmap_data.mat'), 'Z', 'params', '-v7.3');
    fprintf('Saved to: %s\n', out_folder);
end

end

%% =====================================================================
%  LOCAL FUNCTIONS
%  =====================================================================

function plot_4panel_heatmap(vor_pref_smooth, vor_nonpref_smooth, ...
        stat_pref_smooth, stat_nonpref_smooth, ...
        sort_idx, peth_centers, nQual, n_sessions, min_saccades, ...
        FontSize, TitleFontSize, TickFontSize, LineWidth_zero, ...
        doShow, doSavePNG, doSaveFIG, doSaveEPS, out_folder)

    if doShow
        fh = figure('Position', [50 50 1400 1000]);
    else
        fh = figure('Position', [50 50 1400 1000], 'Visible', 'off');
    end
    fh.UserData.date_created        = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fh.UserData.generating_function = 'plot_stationary_saccade_heatmap_ver_1_2';

    % Sort all panels by VOR preferred peak bin
    vor_pref_sorted     = vor_pref_smooth(sort_idx, :);
    vor_nonpref_sorted  = vor_nonpref_smooth(sort_idx, :);
    stat_pref_sorted    = stat_pref_smooth(sort_idx, :);
    stat_nonpref_sorted = stat_nonpref_smooth(sort_idx, :);

    % Common color limits across all 4 panels
    all_vals = [vor_pref_sorted(:); vor_nonpref_sorted(:); ...
                stat_pref_sorted(:); stat_nonpref_sorted(:)];
    clim_range = [0 prctile(all_vals(~isnan(all_vals)), 99)];

    % --- Panel 1: VOR Preferred (top left) ---
    ax1 = subplot(2, 2, 1);
    imagesc(peth_centers, 1:nQual, vor_pref_sorted);
    hold(ax1, 'on');
    yl = ylim(ax1);
    line(ax1, [0 0], yl, 'Color', 'r', 'LineWidth', LineWidth_zero);
    ylabel('NPH unit #', 'FontSize', FontSize);
    title('VOR — Preferred', 'FontSize', TitleFontSize, 'FontWeight', 'bold');
    set(ax1, 'FontSize', TickFontSize, 'XTickLabel', []);
    caxis(ax1, clim_range);

    % --- Panel 2: VOR Non-Preferred (top right) ---
    ax2 = subplot(2, 2, 2);
    imagesc(peth_centers, 1:nQual, vor_nonpref_sorted);
    hold(ax2, 'on');
    yl2 = ylim(ax2);
    line(ax2, [0 0], yl2, 'Color', 'r', 'LineWidth', LineWidth_zero);
    title('VOR — Non-Preferred', 'FontSize', TitleFontSize, 'FontWeight', 'bold');
    set(ax2, 'FontSize', TickFontSize, 'XTickLabel', [], 'YTickLabel', []);
    caxis(ax2, clim_range);

    % --- Panel 3: Spontaneous Preferred (bottom left) ---
    ax3 = subplot(2, 2, 3);
    imagesc(peth_centers, 1:nQual, stat_pref_sorted);
    hold(ax3, 'on');
    yl3 = ylim(ax3);
    line(ax3, [0 0], yl3, 'Color', 'r', 'LineWidth', LineWidth_zero);
    xlabel('Time from saccade (s)', 'FontSize', FontSize);
    ylabel('NPH unit #', 'FontSize', FontSize);
    title('Spontaneous — Preferred', 'FontSize', TitleFontSize, 'FontWeight', 'bold');
    set(ax3, 'FontSize', TickFontSize);
    caxis(ax3, clim_range);

    % --- Panel 4: Spontaneous Non-Preferred (bottom right) ---
    ax4 = subplot(2, 2, 4);
    imagesc(peth_centers, 1:nQual, stat_nonpref_sorted);
    hold(ax4, 'on');
    yl4 = ylim(ax4);
    line(ax4, [0 0], yl4, 'Color', 'r', 'LineWidth', LineWidth_zero);
    xlabel('Time from saccade (s)', 'FontSize', FontSize);
    title('Spontaneous — Non-Preferred', 'FontSize', TitleFontSize, 'FontWeight', 'bold');
    set(ax4, 'FontSize', TickFontSize, 'YTickLabel', []);
    caxis(ax4, clim_range);

    % Shared colorbar
    cb = colorbar(ax4, 'Position', [0.92 0.11 0.02 0.815]);
    cb.Label.String = 'FR (Hz)';
    cb.Label.FontSize = FontSize;
    cb.FontSize = TickFontSize;

    % Supertitle
    sgtitle(sprintf('Saccade PETHs: VOR vs Spontaneous (%d neurons, %d sessions, \\geq%d sacc/dir)', ...
        nQual, n_sessions, min_saccades), 'FontSize', TitleFontSize + 2, 'FontWeight', 'bold');

    % Save
    save_fig_trio(fh, fullfile(out_folder, 'saccade_heatmap_VOR_vs_spont'), ...
        doSavePNG, doSaveFIG, doSaveEPS, doShow);
end

function fr_smooth = compute_peth(myCell, sacc_times, peth_bins, peth_dt, smooth_win, excessBounds)
    nBins_local = length(peth_bins) - 1;
    nEvents = length(sacc_times);
    if nEvents < 1
        fr_smooth = nan(1, nBins_local);
        return
    end
    peth_window_local = [peth_bins(1) peth_bins(end)];
    spk_rel = [];
    for iE = 1:nEvents
        S0 = restrict(myCell, sacc_times(iE)+peth_window_local(1)-excessBounds, ...
                              sacc_times(iE)+peth_window_local(2)+excessBounds);
        S0 = restrict(S0, sacc_times(iE)+peth_window_local(1), ...
                          sacc_times(iE)+peth_window_local(2));
        if ~isempty(S0.t{1})
            spk_rel = [spk_rel; S0.t{1} - sacc_times(iE)]; %#ok<AGROW>
        end
    end
    m_hist = histcounts(spk_rel, peth_bins);
    fr_peth = m_hist / peth_dt / nEvents;
    fr_smooth = smoothdata(movmean(fr_peth, smooth_win), 'gaussian', smooth_win);
end

function bm = compute_burst_mod_value(fr_trace, bl_bins, burst_bins)
    bl_mean = mean(fr_trace(bl_bins), 'omitnan');
    burst_mean = mean(fr_trace(burst_bins), 'omitnan');
    if bl_mean > 0
        bm = (burst_mean - bl_mean) / bl_mean;
    else
        bm = NaN;
    end
end

function out_folder = make_output_folder(save_dir, doSavePNG, doSaveFIG, doSaveEPS)
    run_date_str = datestr(now, 'yyyy-mm-dd');
    run_num = 1;
    out_folder = fullfile(save_dir, sprintf('stat_heatmap_%s_run%d', run_date_str, run_num));
    while exist(out_folder, 'dir')
        run_num = run_num + 1;
        out_folder = fullfile(save_dir, sprintf('stat_heatmap_%s_run%d', run_date_str, run_num));
    end
    if doSavePNG || doSaveFIG || doSaveEPS
        mkdir(out_folder);
    end
end

% =====================================================================
%  VERSION HISTORY
% =====================================================================
% 2026-04-07  v1.2  4-panel layout: VOR Pref/NonPref (top) +
%                    Spont Pref/NonPref (bottom). All panels sorted by
%                    VOR preferred peak bin. Removed local save_fig_trio
%                    and find_cell_index (use 06_utilities/ versions).
%                    Shared color limits across all 4 panels (99th pctile).
%                    Versioned filename. min_saccades default = 10.
%                    Workspace mode re-thresholding stores VOR data.
% 2026-04-04  v1.1  Log-scale option, workspace re-thresholding, VOR
%                    data stored in Z for reference.
% 2026-04-04  v1.0  Initial: 2-panel (stationary pref/nonpref only).
