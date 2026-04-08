function [Z, params] = plot_stationary_saccade_heatmap(varargin)
% plot_stationary_saccade_heatmap  v1.1
% Population heatmap of saccade PETHs during stationary (spontaneous) periods.
% Analog of Figure 3C but restricted to platform-stationary saccades.
% Log-scaled colormap with Hz tick labels. Both panels sorted by
% preferred peak bin (matching row order).
%
%   Uses Chenhao's saccade onset times (saccade_ONSET_all_unit.mat).
%   Preferred direction determined from VOR saccade burst modulation.
%   Two panels: preferred | non-preferred, rows sorted by preferred peak bin.
%
%   INPUTS (name-value pairs via varargin):
%     tfilelist_path    - path to NPH_all_tfilelist.mat (default: standard)
%     min_saccades      - min stationary saccades per dir (default: 5)
%     peth_window       - [start end] in s (default: [-0.2 0.2])
%     peth_dt           - bin width in s (default: 0.005)
%     peth_smooth_win   - smoothing window in bins (default: 5)
%     burst_window      - for preferred dir detection (default: [-0.05 0.05])
%     baseline_window   - for preferred dir detection (default: [-0.2 -0.1])
%     doSavePNG/FIG/EPS - save flags (default: 1)
%     doShow            - display figure (default: 1)
%     save_dir          - output root (default: ___analysis_output)
%
%   OUTPUTS:
%     Z      - struct with FR matrices, sort indices, bin centers (like Jeff's Z)
%     params - struct of all free parameter values
%
%   DEPENDENCIES:
%     LoadSpikesQuick, SelectTS, restrict, pushdir, popdir, HD_GetSSN,
%     process_varargin, load_chenhao_saccades

%% =====================================================================
%  DEFAULTS
%  =====================================================================
tfilelist_path    = 'C:\Jeff\U01\datatouse\NPH all tfilelist.mat';
save_dir          = 'C:\Jeff\U01\datatouse\___analysis_output';
precomputed_Z     = [];   % pass Z struct from prior run to skip recomputation

% PETH parameters (matching v1.3 / Jeff's heatmap code)
peth_window       = [-0.2 0.2];
peth_dt           = 0.005;
peth_smooth_win   = 5;
peth_excessBounds = 1;

% Burst modulation windows (for preferred direction detection)
burst_window      = [-0.05 0.05];
baseline_window   = [-0.2 -0.1];

% Thresholds
min_saccades      = 5;

% Figure
FontSize          = 28;
TitleFontSize     = 24;
TickFontSize      = 18;
LineWidth_zero    = 2;
doLogScale        = 0;   % 0 = linear (default), 1 = log10 colormap

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
params.version          = '1.1';
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
if ~isempty(precomputed_Z) && isfield(precomputed_Z, 'FR_preferred_smooth')
    fprintf('=== WORKSPACE MODE (Z struct with %d neurons) ===\n', precomputed_Z.n_neurons);
    Z = precomputed_Z;

    % Allow re-thresholding from workspace
    if isfield(Z, 'all_FR_stat_pref') && min_saccades ~= Z.min_saccades
        fprintf('  Re-thresholding: %d -> %d min saccades\n', Z.min_saccades, min_saccades);
        new_keep = Z.all_n_stat_CCW >= min_saccades & Z.all_n_stat_CW >= min_saccades & Z.all_qualifies;
        Z.FR_preferred_smooth   = smoothdata(Z.all_FR_stat_pref(new_keep, :), 2);
        Z.FR_nonpreferred_smooth = smoothdata(Z.all_FR_stat_nonpref(new_keep, :), 2);
        Z.FR_preferred          = Z.all_FR_stat_pref(new_keep, :);
        Z.FR_nonpreferred       = Z.all_FR_stat_nonpref(new_keep, :);
        Z.neuron_names          = Z.all_neuron_names(new_keep);
        Z.n_neurons             = sum(new_keep);
        Z.min_saccades          = min_saccades;
        [~, peak_bin_re] = max(Z.FR_preferred, [], 2);
        [~, Z.sPreferred] = sort(peak_bin_re);
    end

    % Output folder
    run_date_str = datestr(now, 'yyyy-mm-dd');
    run_num = 1;
    out_folder = fullfile(save_dir, sprintf('stat_heatmap_%s_run%d', run_date_str, run_num));
    while exist(out_folder, 'dir')
        run_num = run_num + 1;
        out_folder = fullfile(save_dir, sprintf('stat_heatmap_%s_run%d', run_date_str, run_num));
    end
    if doSavePNG || doSaveFIG || doSaveEPS; mkdir(out_folder); end
    params.output_folder = out_folder;

    % Jump to plotting
    nQual = Z.n_neurons;
    n_sessions = Z.n_sessions;
    peth_centers = Z.binCenters;
    qual_stat_pref_smooth    = Z.FR_preferred_smooth;
    qual_stat_nonpref_smooth = Z.FR_nonpreferred_smooth;
    sort_pref = Z.sPreferred;

    % --- PLOT (workspace mode uses same code as full run) ---
    plot_heatmap_figure(qual_stat_pref_smooth, qual_stat_nonpref_smooth, sort_pref, ...
        peth_centers, nQual, n_sessions, min_saccades, ...
        FontSize, TitleFontSize, TickFontSize, LineWidth_zero, doLogScale, doShow, ...
        doSavePNG, doSaveFIG, doSaveEPS, out_folder);

    if doSavePNG || doSaveFIG || doSaveEPS
        save(fullfile(out_folder, 'heatmap_data.mat'), 'Z', 'params');
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
run_date_str = datestr(now, 'yyyy-mm-dd');
run_num = 1;
out_folder = fullfile(save_dir, sprintf('stat_heatmap_%s_run%d', run_date_str, run_num));
while exist(out_folder, 'dir')
    run_num = run_num + 1;
    out_folder = fullfile(save_dir, sprintf('stat_heatmap_%s_run%d', run_date_str, run_num));
end
if doSavePNG || doSaveFIG || doSaveEPS; mkdir(out_folder); end
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
    [sess_path, neuron_name] = fileparts(tfile_path);
    neuron_names{iN} = neuron_name;

    if verbose && mod(iN, 10) == 0
        fprintf('  %d/%d ...\n', iN, nNeurons);
    end

    try
        pushdir(sess_path);
        SSN = HD_GetSSN;
        session_names{iN} = SSN;

        %% Load spikes (only reload if new session)
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

            % Split VOR vs stationary
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
        cell_idx = find_cell_index(spike_cell_names, neuron_name, SSN);
        if isempty(cell_idx); popdir; continue; end
        myCell = SelectTS([], S_all, cell_idx);

        n_stat_CCW(iN) = length(sess_temporal_stat);
        n_stat_CW(iN)  = length(sess_nasal_stat);
        n_vor_CCW(iN)  = length(sess_temporal_move);
        n_vor_CW(iN)   = length(sess_nasal_move);

        %% Compute 4 PETHs
        fr_vor_CCW  = compute_peth(myCell, sess_temporal_move, peth_bins, peth_dt, peth_smooth_win, peth_excessBounds);
        fr_vor_CW   = compute_peth(myCell, sess_nasal_move, peth_bins, peth_dt, peth_smooth_win, peth_excessBounds);
        fr_stat_CCW = compute_peth(myCell, sess_temporal_stat, peth_bins, peth_dt, peth_smooth_win, peth_excessBounds);
        fr_stat_CW  = compute_peth(myCell, sess_nasal_stat, peth_bins, peth_dt, peth_smooth_win, peth_excessBounds);

        %% Determine preferred direction from VOR burst modulation
        bm_CCW = compute_burst_mod_trace(fr_vor_CCW, base_mask, burst_mask);
        bm_CW  = compute_burst_mod_trace(fr_vor_CW, base_mask, burst_mask);

        if isnan(bm_CCW) && isnan(bm_CW)
            pref_dir = 'none';
        elseif isnan(bm_CW) || bm_CCW >= bm_CW
            pref_dir = 'CCW';
        else
            pref_dir = 'CW';
        end
        pref_dirs{iN} = pref_dir;

        if strcmp(pref_dir, 'CCW')
            FR_stat_pref(iN,:)    = fr_stat_CCW;
            FR_stat_nonpref(iN,:) = fr_stat_CW;
            FR_vor_pref(iN,:)     = fr_vor_CCW;
            FR_vor_nonpref(iN,:)  = fr_vor_CW;
        else
            FR_stat_pref(iN,:)    = fr_stat_CW;
            FR_stat_nonpref(iN,:) = fr_stat_CCW;
            FR_vor_pref(iN,:)     = fr_vor_CW;
            FR_vor_nonpref(iN,:)  = fr_vor_CCW;
        end

        %% Check if qualifies
        qualifies(iN) = n_stat_CCW(iN) >= min_saccades && n_stat_CW(iN) >= min_saccades ...
            && ~strcmp(pref_dir, 'none');

        popdir;
    catch ME
        if verbose; fprintf('  ERROR %s: %s\n', neuron_name, ME.message); end
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

%% =====================================================================
%  PROCESS (matching Jeff's heatmap code)
%  =====================================================================
% Smooth across bins
qual_stat_pref_smooth    = smoothdata(qual_stat_pref, 2);
qual_stat_nonpref_smooth = smoothdata(qual_stat_nonpref, 2);

% Sort by peak bin of preferred
[~, peak_bin_pref] = max(qual_stat_pref, [], 2);
[~, sort_pref]     = sort(peak_bin_pref);

%% =====================================================================
%  BUILD Z STRUCT (matching Jeff's convention)
%  =====================================================================
Z = struct();
Z.FR_preferred          = qual_stat_pref;
Z.FR_nonpreferred       = qual_stat_nonpref;
Z.FR_preferred_smooth   = qual_stat_pref_smooth;
Z.FR_nonpreferred_smooth = qual_stat_nonpref_smooth;
Z.sPreferred            = sort_pref;
Z.binCenters            = peth_centers;
Z.peth_bins             = peth_bins;
Z.neuron_names          = neuron_names(keep_idx);
Z.session_names         = session_names(keep_idx);
Z.pref_dirs             = pref_dirs(keep_idx);
Z.n_stat_CCW            = n_stat_CCW(keep_idx);
Z.n_stat_CW             = n_stat_CW(keep_idx);
Z.n_neurons             = nQual;
Z.n_sessions            = n_sessions;
Z.min_saccades          = min_saccades;

% Also store VOR data for reference
Z.FR_vor_pref           = FR_vor_pref(keep_idx, :);
Z.FR_vor_nonpref        = FR_vor_nonpref(keep_idx, :);

% All neuron data (unfiltered) for re-thresholding
Z.all_FR_stat_pref      = FR_stat_pref;
Z.all_FR_stat_nonpref   = FR_stat_nonpref;
Z.all_qualifies         = qualifies;
Z.all_neuron_names      = neuron_names;
Z.all_n_stat_CCW        = n_stat_CCW;
Z.all_n_stat_CW         = n_stat_CW;

%% =====================================================================
%  PLOT + SAVE
%  =====================================================================
plot_heatmap_figure(qual_stat_pref_smooth, qual_stat_nonpref_smooth, sort_pref, ...
    peth_centers, nQual, n_sessions, min_saccades, ...
    FontSize, TitleFontSize, TickFontSize, LineWidth_zero, doLogScale, doShow, ...
    doSavePNG, doSaveFIG, doSaveEPS, out_folder);

if doSavePNG || doSaveFIG || doSaveEPS
    save(fullfile(out_folder, 'heatmap_data.mat'), 'Z', 'params');
    fprintf('Saved to: %s\n', out_folder);
end

end

%% =====================================================================
%  LOCAL HELPERS
%  =====================================================================

function plot_heatmap_figure(pref_smooth, nonpref_smooth, sort_idx, ...
        peth_centers, nQual, n_sessions, min_saccades, ...
        FontSize, TitleFontSize, TickFontSize, LineWidth_zero, doLogScale, doShow, ...
        doSavePNG, doSaveFIG, doSaveEPS, out_folder)

    if doShow; fh = figure('Position', [50 50 1200 700]);
    else; fh = figure('Position', [50 50 1200 700], 'Visible', 'off'); end
    fh.UserData.date_created        = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fh.UserData.generating_function = 'plot_stationary_saccade_heatmap v1.2';

    % Prepare data (sorted by preferred peak bin for BOTH panels)
    pref_sorted    = pref_smooth(sort_idx, :);
    nonpref_sorted = nonpref_smooth(sort_idx, :);

    if doLogScale
        log_eps = 1;
        plot_pref    = log10(pref_sorted + log_eps);
        plot_nonpref = log10(nonpref_sorted + log_eps);
        clim_range = [0 max([plot_pref(:); plot_nonpref(:)])];
        cb_hz_ticks = [0 10 30 50 100 200 300];
        cb_hz_ticks = cb_hz_ticks(cb_hz_ticks <= 10^clim_range(2));
        cb_log_ticks = log10(cb_hz_ticks + log_eps);
    else
        plot_pref    = pref_sorted;
        plot_nonpref = nonpref_sorted;
        clim_range = [0 max([plot_pref(:); plot_nonpref(:)])];
    end

    % Left panel: preferred
    ax1 = subplot(1, 2, 1);
    imagesc(peth_centers, 1:nQual, plot_pref);
    hold(ax1, 'on');
    al = axis(ax1);
    line(ax1, [0 0], [al(3) al(4)], 'Color', 'r', 'LineWidth', LineWidth_zero);
    xlabel('Time (s)', 'FontSize', FontSize)
    ylabel('NPH unit #', 'FontSize', FontSize)
    title('Preferred Saccade Direction', 'FontSize', TitleFontSize, 'FontWeight', 'bold')
    set(ax1, 'FontSize', TickFontSize);
    caxis(ax1, clim_range);

    % Right panel: non-preferred (SAME row order)
    ax2 = subplot(1, 2, 2);
    imagesc(peth_centers, 1:nQual, plot_nonpref);
    hold(ax2, 'on');
    al2 = axis(ax2);
    line(ax2, [0 0], [al2(3) al2(4)], 'Color', 'r', 'LineWidth', LineWidth_zero);
    xlabel('Time (s)', 'FontSize', FontSize)
    title('Non-Preferred Saccade Direction', 'FontSize', TitleFontSize, 'FontWeight', 'bold')
    set(ax2, 'FontSize', TickFontSize);
    caxis(ax2, clim_range);

    % Colorbar
    cb = colorbar(ax2);
    cb.Label.FontSize = FontSize;
    if doLogScale
        cb.Label.String = 'FR (Hz)';
        cb.Ticks = cb_log_ticks;
        cb.TickLabels = arrayfun(@(x) sprintf('%d', x), cb_hz_ticks, 'UniformOutput', false);
    else
        cb.Label.String = 'FR (Hz)';
    end

    sgtitle(sprintf('Spontaneous Saccade PETHs (%d neurons, %d sessions, >= %d saccades/dir)', ...
        nQual, n_sessions, min_saccades), 'FontSize', TitleFontSize + 2, 'FontWeight', 'bold');

    save_fig_trio(fh, fullfile(out_folder, 'stationary_saccade_heatmap'), ...
        doSavePNG, doSaveFIG, doSaveEPS, doShow);
end

function fr_smooth = compute_peth(myCell, sacc_times, peth_bins, peth_dt, smooth_win, excessBounds)
    nBins = length(peth_bins) - 1;
    nEvents = length(sacc_times);
    if nEvents < 1
        fr_smooth = nan(1, nBins);
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

function bm = compute_burst_mod_trace(fr_trace, bl_bins, burst_bins)
    bl_mean = mean(fr_trace(bl_bins), 'omitnan');
    burst_mean = mean(fr_trace(burst_bins), 'omitnan');
    if bl_mean > 0; bm = (burst_mean - bl_mean) / bl_mean; else; bm = NaN; end
end

function cell_idx = find_cell_index(spike_cell_names, neuron_name, SSN)
    cell_idx = [];
    for cc = 1:length(spike_cell_names)
        if strcmp(spike_cell_names{cc}, neuron_name); cell_idx = cc; return; end
    end
    for cc = 1:length(spike_cell_names)
        if contains(spike_cell_names{cc}, strrep(neuron_name, [SSN '-'], ''))
            cell_idx = cc; return;
        end
    end
end

function save_fig_trio(fh, filepath_base, doSavePNG, doSaveFIG, doSaveEPS, doShow)
    if doSavePNG; print(fh, [filepath_base '.png'], '-dpng', '-r300'); end
    if doSaveFIG; savefig(fh, [filepath_base '.fig']); end
    if doSaveEPS; print(fh, [filepath_base '.eps'], '-depsc2'); end
    if ~doShow; close(fh); end
end
