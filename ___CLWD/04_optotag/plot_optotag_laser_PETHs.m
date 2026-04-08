function [fig_handles, params, cell_data] = plot_optotag_laser_PETHs(cfg_in, varargin)
% plot_optotag_laser_PETHs  Figure S4: Laser PETHs for all 7 optotagged neurons
%
%   Produces two figure versions:
%     Version A: 8-panel layout (7 laser PETHs + AHV heatmap in slot 8)
%     Version B: 7-panel layout with AHV TC inset in each panel
%
%   Each PETH panel shows:
%     - Top: spike raster aligned to laser onset (black dots)
%     - Bottom: mean firing rate PETH (red line)
%     - Blue shading for laser-on interval
%     - Green line during laser interval (placeholder for significance)
%
%   If cfg.cell_data is provided, skips data loading and replots from cache.
%
%   INPUTS:
%       cfg_in      - config structure (all fields optional)
%       varargin    - name-value pairs overriding any cfg field
%
%   CFG FIELDS:
%       .cell_data      - precomputed cell data (skip reloading if provided)
%       .window         - PETH window in s (default: [-1 2])
%       .dt             - FR bin size in s (default: 0.05)
%       .laser_dur      - laser pulse duration in s (default: 1.0)
%       .smoothFR       - smooth FR trace (default: 1)
%       .smooth_win     - Gaussian smoothing window in bins (default: 3)
%       .AHV_binEdges   - AHV TC bin edges (default: linspace(-200,200,67))
%       .AHV_minOcc     - AHV TC min occupancy in s (default: 0.5)
%       .smoothTC_win   - AHV TC smoothing window (default: 5)
%       .FontSize       - base font size (default: 24)
%       .rasterTickSize - raster dot size (default: 3)
%       .showTitles     - show cell name titles (default: 0)
%       .figWidth       - figure width in pixels (default: 2200)
%       .figHeight      - figure height in pixels (default: 1200)
%       .doShow         - display figures (default: 0)
%       .doSavePNG      - save .png (default: 1)
%       .doSaveFIG      - save .fig (default: 0)
%       .savePath       - root save path
%       .shuff_results  - cell array from compute_optotag_shuffle (optional;
%                         when provided, draws real significance bars and
%                         shuffle mean +/- SEM instead of placeholder)
%       .verbose        - print progress (default: 1)
%
%   OUTPUTS:
%       fig_handles - [fig_A, fig_B] figure handles
%       params      - parameter struct
%       cell_data   - cell array of per-neuron data structs
%
%   DEPENDENCIES:
%       LoadSpikesQuick, SelectTS, SortBrainstemEventLabels3, HD_GetSSN,
%       get_AHVtsd, TuningCurves, tsd, LoadSpikes, LoadEvents,
%       ProcessConfig2, process_varargin, pushdir, popdir
%

%% Config defaults
cfg_def = [];
cfg_def.cell_data       = [];
cfg_def.window          = [-1 2];
cfg_def.dt              = 0.05;
cfg_def.laser_dur       = 1.0;
cfg_def.smoothFR        = 1;
cfg_def.smooth_win      = 3;
cfg_def.AHV_binEdges    = linspace(-200, 200, 67);  % 6 deg/s bins
cfg_def.AHV_minOcc      = 0.5;
cfg_def.smoothTC_win    = 5;
cfg_def.FontSize        = 24;
cfg_def.rasterTickSize  = 3;
cfg_def.showTitles      = 0;
cfg_def.figWidth        = 2200;
cfg_def.figHeight       = 1200;
cfg_def.doShow          = 0;
cfg_def.doSavePNG       = 1;
cfg_def.doSaveFIG       = 0;
cfg_def.doVersionA      = 0;     % PETHs only (no TCs, no heatmap)
cfg_def.doVersionB      = 0;     % TC insets + heatmap
cfg_def.doVersionC      = 1;     % PETHs + heatmap (default for paper)
cfg_def.savePath        = 'C:\Jeff\U01\datatouse\___analysis_output';
cfg_def.shuff_results   = [];    % cell array from compute_optotag_shuffle (one per cell)
cfg_def.verbose         = 1;

if nargin < 1 || isempty(cfg_in); cfg_in = []; end
cfg = ProcessConfig2(cfg_def, cfg_in);
process_varargin(varargin);

%% Store params
params                  = cfg;
params.function_name    = 'plot_optotag_laser_PETHs';
params.version          = '2.4';
params.date_run         = datestr(now, 'yyyy-mm-dd');

%% Define the 7 optotagged neurons
data_root = 'C:\Jeff\U01\datatouse';
opto_cells = {
    fullfile(data_root, 'M271', 'M271-2021-08-30'),    'TT02_1',  'M271-2021-08-30-TT02-1';
    fullfile(data_root, 'M282', 'M282-2022-02-10'),     'TT03_2',  'M282-2022-02-10-TT03-2';
    fullfile(data_root, 'M282', 'M282-2022-02-15'),     'TT07_1',  'M282-2022-02-15-TT07-1';
    fullfile(data_root, 'M391', 'M391-2022-12-07-1'),   'TT01_1',  'M391-2022-12-07-1-TT01-1';
    fullfile(data_root, 'M391', 'M391-2022-12-07-2'),   'TT02_1',  'M391-2022-12-07-2-TT02-1';
    fullfile(data_root, 'M402', 'M402-2022-11-02-2'),   'TT02_3',  'M402-2022-11-02-2-TT02-3';
    fullfile(data_root, 'M444', 'M444-2023-08-18'),     'TT07_1',  'M444-2023-08-18-TT07-1';
};
nCells = size(opto_cells, 1);

%% ========================================================================
%  DATA COMPUTATION (skip if cell_data provided)
%  ========================================================================
if ~isempty(cfg.cell_data)
    cell_data = cfg.cell_data;
    if cfg.verbose; fprintf('Using precomputed cell_data (%d cells).\n', length(cell_data)); end
else
    cell_data = cell(nCells, 1);

    for iC = 1:nCells
        sess_dir    = opto_cells{iC, 1};
        tfile_tag   = opto_cells{iC, 2};
        cell_label  = opto_cells{iC, 3};

        if cfg.verbose
            fprintf('\n=== Cell %d/%d: %s ===\n', iC, nCells, cell_label);
        end

        pushdir(sess_dir);
        SSN = HD_GetSSN;

        %% Load spikes
        [S_all, cell_names_all, spk_params] = LoadSpikesQuick();

        %% Get session duration (for circular shuffle)
        events_ts = LoadEvents([]);
        stop_label_idx = [];
        for iL = 1:length(events_ts.label)
            if strcmp(events_ts.label{iL}, 'Stopping Recording')
                stop_label_idx = iL; break
            end
        end
        if isempty(stop_label_idx)
            warning('Could not find Stopping Recording event for %s.', SSN);
            sess_duration = max(cellfun(@max, S_all.t));  % fallback: last spike
        else
            sess_duration = events_ts.t{stop_label_idx}(end) - spk_params.starttime;
        end

        %% Get laser times
        [~, ~, laser_on, ~] = SortBrainstemEventLabels3;
        if isnan(laser_on)
            warning('No laser events for %s. Skipping.', SSN);
            popdir; continue
        end
        nTrials = length(laser_on);
        if cfg.verbose; fprintf('  %d laser pulses found.\n', nTrials); end

        %% Match cell
        cell_idx = [];
        for kk = 1:length(cell_names_all)
            if contains(cell_names_all{kk}, tfile_tag)
                cell_idx = kk; break
            end
        end
        if isempty(cell_idx)
            warning('Could not match %s. Skipping.', tfile_tag);
            popdir; continue
        end

        myCell = SelectTS([], S_all, cell_idx);
        spike_times = myCell.t{1};

        %% Build spike raster and FR PETH
        bin_edges   = cfg.window(1) : cfg.dt : cfg.window(2);
        bin_centers = bin_edges(1:end-1) + cfg.dt / 2;
        nBins       = length(bin_centers);

        raster_spk_times = [];
        raster_trial_num = [];
        fr_hist          = zeros(1, nBins);

        for iT = 1:nTrials
            t0 = laser_on(iT);
            spk_in_win = spike_times(spike_times >= t0 + cfg.window(1) & ...
                                      spike_times <  t0 + cfg.window(2));
            spk_rel = spk_in_win - t0;
            raster_spk_times = [raster_spk_times; spk_rel(:)];
            raster_trial_num = [raster_trial_num; repmat(iT, length(spk_rel), 1)];
            if ~isempty(spk_rel)
                fr_hist = fr_hist + histcounts(spk_rel, bin_edges);
            end
        end

        fr_mean_vals = fr_hist / (nTrials * cfg.dt);
        if cfg.smoothFR && cfg.smooth_win > 1
            fr_smooth_vals = smoothdata(fr_mean_vals, 'gaussian', cfg.smooth_win);
        else
            fr_smooth_vals = fr_mean_vals;
        end

        %% Compute AHV tuning curve
        if cfg.verbose; fprintf('  Computing AHV tuning curve...\n'); end
        AHV_tsd = get_AHVtsd([]);

        AHV_dt = median(diff(AHV_tsd.tvec));
        cfg_tc          = [];
        cfg_tc.nBins    = length(cfg.AHV_binEdges) - 1;
        cfg_tc.binEdges = {cfg.AHV_binEdges};
        cfg_tc.minOcc   = cfg.AHV_minOcc;
        cfg_tc.occ_dt   = AHV_dt;
        tc_AHV          = TuningCurves(cfg_tc, myCell, AHV_tsd);

        tc_vals = tc_AHV.tc(1,:);
        tc_occ  = tc_AHV.occ_hist;
        AHV_bin_centers = cfg.AHV_binEdges(1:end-1) + diff(cfg.AHV_binEdges(1:2))/2;
        good_bins_ahv = tc_occ > cfg.AHV_minOcc;
        tc_smooth_ahv = nan(size(tc_vals));
        if sum(good_bins_ahv) > 2
            tc_smooth_ahv(good_bins_ahv) = smoothdata(tc_vals(good_bins_ahv), 'gaussian', cfg.smoothTC_win);
        end

        %% Store
        cell_data{iC} = struct();
        cell_data{iC}.cell_label       = cell_label;
        cell_data{iC}.SSN              = SSN;
        cell_data{iC}.tfile_tag        = tfile_tag;
        cell_data{iC}.nTrials          = nTrials;
        cell_data{iC}.laser_on         = laser_on;
        cell_data{iC}.raster_spk_times = raster_spk_times;
        cell_data{iC}.raster_trial_num = raster_trial_num;
        cell_data{iC}.bin_centers      = bin_centers;
        cell_data{iC}.fr_mean          = fr_mean_vals;
        cell_data{iC}.fr_smooth        = fr_smooth_vals;
        cell_data{iC}.AHV_bin_centers  = AHV_bin_centers;
        cell_data{iC}.AHV_tc_raw       = tc_vals;
        cell_data{iC}.AHV_tc_smooth    = tc_smooth_ahv;
        cell_data{iC}.AHV_occ          = tc_occ;
        cell_data{iC}.AHV_good_bins    = good_bins_ahv;
        cell_data{iC}.spike_times      = spike_times;
        cell_data{iC}.session_duration = sess_duration;

        popdir;
    end
end

%% ========================================================================
%  CREATE SAVE FOLDER (new folder per run)
%  ========================================================================
base_folder_name = ['figS4_optotag_PETHs_' params.date_run];
run_num = 1;
save_subfolder = fullfile(cfg.savePath, [base_folder_name '_run' num2str(run_num)]);
while exist(save_subfolder, 'dir')
    run_num = run_num + 1;
    save_subfolder = fullfile(cfg.savePath, [base_folder_name '_run' num2str(run_num)]);
end
mkdir(save_subfolder);
params.save_subfolder = save_subfolder;
if cfg.verbose; fprintf('Save folder: %s\n', save_subfolder); end

% Always save .mat
fname_mat = fullfile(save_subfolder, 'figS4_optotag_data.mat');
save(fname_mat, 'cell_data', 'params', '-v7.3');
if cfg.verbose; fprintf('\nSaved MAT: %s\n', fname_mat); end

%% ========================================================================
%  FIGURE VERSION A: 7 laser PETHs only (no TCs, no heatmap)
%  ========================================================================
fig_A = [];
if cfg.doVersionA
fig_A = make_figure_A(cfg, cell_data, nCells, params);
fname_png_A = fullfile(save_subfolder, 'figS4_versionA_PETHs_only.png');
set(fig_A, 'Visible', 'on'); drawnow;
print(fig_A, fname_png_A, '-dpng', '-r300');
if cfg.verbose; fprintf('Saved PNG: %s\n', fname_png_A); end
if cfg.doSaveFIG
    savefig(fig_A, fullfile(save_subfolder, 'figS4_versionA_PETHs_only.fig'));
end
if ~cfg.doShow; close(fig_A); end
end

%% ========================================================================
%  FIGURE VERSION B: 7 panels with AHV TC inset + heatmap in slot 8
%  ========================================================================
fig_B = [];
if cfg.doVersionB
fig_B = make_figure_B(cfg, cell_data, nCells, params);
fname_png_B = fullfile(save_subfolder, 'figS4_versionB_TC_inset.png');
set(fig_B, 'Visible', 'on'); drawnow;
print(fig_B, fname_png_B, '-dpng', '-r300');
if cfg.verbose; fprintf('Saved PNG: %s\n', fname_png_B); end
if cfg.doSaveFIG
    savefig(fig_B, fullfile(save_subfolder, 'figS4_versionB_TC_inset.fig'));
end
if ~cfg.doShow; close(fig_B); end
end

%% ========================================================================
%  FIGURE VERSION C: 7 laser PETHs + heatmap in slot 8 (no TC insets)
%  ========================================================================
fig_C = [];
if cfg.doVersionC
fig_C = make_figure_C(cfg, cell_data, nCells, params);
fname_png_C = fullfile(save_subfolder, 'figS4_versionC_heatmap.png');
set(fig_C, 'Visible', 'on'); drawnow;
print(fig_C, fname_png_C, '-dpng', '-r300');
if cfg.verbose; fprintf('Saved PNG: %s\n', fname_png_C); end
if cfg.doSaveFIG
    savefig(fig_C, fullfile(save_subfolder, 'figS4_versionC_heatmap.fig'));
end
if ~cfg.doShow; close(fig_C); end
end

fig_handles = [fig_A, fig_B, fig_C];
if cfg.verbose; fprintf('\nDone.\n'); end
end

%% ========================================================================
%  LOCAL FUNCTION: FIGURE VERSION A (PETHs only, no TCs, no heatmap)
%  ========================================================================
function fig_h = make_figure_A(cfg, cell_data, nCells, params)

fig_h = figure('Position', [50 50 cfg.figWidth cfg.figHeight], ...
    'Color', 'w', 'Visible', 'off');
fig_h.UserData.function_name = params.function_name;
fig_h.UserData.version       = params.version;
fig_h.UserData.date_created  = params.date_run;
fig_h.UserData.fig_type      = 'versionA_PETHs_only';

nCols = 4; nRows = 2;
panel_col = [1 2 3 4  1 2 3];
panel_row = [1 1 1 1  2 2 2];

margin_left = 0.06; margin_right = 0.02;
margin_top = 0.04; margin_bottom = 0.08;
hspace = 0.04; vspace = 0.08;

panel_w = (1 - margin_left - margin_right - (nCols-1)*hspace) / nCols;
panel_h = (1 - margin_top - margin_bottom - (nRows-1)*vspace) / nRows;
raster_frac = 0.55; fr_frac = 0.40; gap_frac = 0.05;

for iC = 1:nCells
    if isempty(cell_data{iC}); continue; end
    cd_cell = cell_data{iC};
    col = panel_col(iC); row = panel_row(iC);
    x0 = margin_left + (col-1)*(panel_w + hspace);
    y0_panel = margin_bottom + (nRows - row)*(panel_h + vspace);

    raster_h = panel_h * raster_frac;
    fr_h     = panel_h * fr_frac;
    gap_h    = panel_h * gap_frac;

    ax_r = axes('Position', [x0, y0_panel + fr_h + gap_h, panel_w, raster_h]);
    ax_f = axes('Position', [x0, y0_panel, panel_w, fr_h]);

    plot_raster_panel(ax_r, cd_cell, cfg, iC);
    plot_fr_panel(ax_f, cd_cell, cfg, iC);
end

end

%% ========================================================================
%  LOCAL FUNCTION: FIGURE VERSION B (AHV TC inset per panel + heatmap)
%  ========================================================================
function fig_h = make_figure_B(cfg, cell_data, nCells, params)

fig_h = figure('Position', [50 50 cfg.figWidth cfg.figHeight+200], ...
    'Color', 'w', 'Visible', 'off');
fig_h.UserData.function_name = params.function_name;
fig_h.UserData.version       = params.version;
fig_h.UserData.date_created  = params.date_run;
fig_h.UserData.fig_type      = 'versionB_TC_inset';

nCols = 4; nRows = 2;
panel_col = [1 2 3 4  1 2 3];
panel_row = [1 1 1 1  2 2 2];

margin_left = 0.06; margin_right = 0.02;
margin_top = 0.03; margin_bottom = 0.06;
hspace = 0.04; vspace = 0.06;

panel_w = (1 - margin_left - margin_right - (nCols-1)*hspace) / nCols;
panel_h = (1 - margin_top - margin_bottom - (nRows-1)*vspace) / nRows;
raster_frac = 0.38; fr_frac = 0.30; tc_frac = 0.25; gap_frac = 0.035;

for iC = 1:nCells
    if isempty(cell_data{iC}); continue; end
    cd_cell = cell_data{iC};
    col = panel_col(iC); row = panel_row(iC);
    x0 = margin_left + (col-1)*(panel_w + hspace);
    y0_panel = margin_bottom + (nRows - row)*(panel_h + vspace);

    raster_h = panel_h * raster_frac;
    fr_h     = panel_h * fr_frac;
    tc_h     = panel_h * tc_frac;
    gap_h    = panel_h * gap_frac;

    y_tc     = y0_panel;
    y_fr     = y0_panel + tc_h + gap_h;
    y_raster = y_fr + fr_h + gap_h;

    ax_r  = axes('Position', [x0, y_raster, panel_w, raster_h]);
    ax_f  = axes('Position', [x0, y_fr, panel_w, fr_h]);
    ax_tc = axes('Position', [x0, y_tc, panel_w, tc_h]);

    plot_raster_panel(ax_r, cd_cell, cfg, iC);
    plot_fr_panel(ax_f, cd_cell, cfg, iC);

    % In Version B, remove the xlabel from FR panel (TC panel below has its own)
    xlabel(ax_f, '');
    set(ax_f, 'XTickLabel', []);

    % AHV TC panel
    axes(ax_tc); hold on;
    good = cd_cell.AHV_good_bins;
    if any(good)
        plot(cd_cell.AHV_bin_centers(good), cd_cell.AHV_tc_smooth(good), 'k', 'LineWidth', 2);
    end
    tc_max = max(cd_cell.AHV_tc_smooth(good));
    if isempty(tc_max) || tc_max == 0; tc_max = 10; end
    set(ax_tc, 'XLim', [-200 200], 'YLim', [0 tc_max*1.2], ...
        'TickDir', 'out', 'FontSize', cfg.FontSize - 6, ...
        'FontName', 'Arial', ...
        'XTick', [-200 0 200], 'Box', 'off');
    if iC == 1 || iC == 5
        ylabel('FR', 'FontSize', cfg.FontSize - 6, 'FontName', 'Arial');
    end
    if iC == 5
        xlabel('AHV (deg/s)', 'FontSize', cfg.FontSize - 6);
    end

    % Add x-axis label to FR panel for panel 5
    if iC == 5
        set(ax_f, 'XTickLabel', {'-1','0','1','2'});
    end
end

% Panel 8: AHV heatmap
col = 4; row = 2;
x0 = margin_left + (col-1)*(panel_w + hspace);
y0_panel = margin_bottom + (nRows - row)*(panel_h + vspace);
ax_heat = axes('Position', [x0, y0_panel, panel_w, panel_h]);
plot_ahv_heatmap(ax_heat, cell_data, cfg);

end

%% ========================================================================
%  LOCAL FUNCTION: FIGURE VERSION C (PETHs + heatmap, no TC insets)
%  ========================================================================
function fig_h = make_figure_C(cfg, cell_data, nCells, params)

fig_h = figure('Position', [50 50 cfg.figWidth cfg.figHeight], ...
    'Color', 'w', 'Visible', 'off');
fig_h.UserData.function_name = params.function_name;
fig_h.UserData.version       = params.version;
fig_h.UserData.date_created  = params.date_run;
fig_h.UserData.fig_type      = 'versionC_heatmap';

nCols = 4; nRows = 2;
panel_col = [1 2 3 4  1 2 3];
panel_row = [1 1 1 1  2 2 2];

margin_left = 0.06; margin_right = 0.02;
margin_top = 0.04; margin_bottom = 0.08;
hspace = 0.04; vspace = 0.08;

panel_w = (1 - margin_left - margin_right - (nCols-1)*hspace) / nCols;
panel_h = (1 - margin_top - margin_bottom - (nRows-1)*vspace) / nRows;
raster_frac = 0.55; fr_frac = 0.40; gap_frac = 0.05;

for iC = 1:nCells
    if isempty(cell_data{iC}); continue; end
    cd_cell = cell_data{iC};
    col = panel_col(iC); row = panel_row(iC);
    x0 = margin_left + (col-1)*(panel_w + hspace);
    y0_panel = margin_bottom + (nRows - row)*(panel_h + vspace);

    raster_h = panel_h * raster_frac;
    fr_h     = panel_h * fr_frac;
    gap_h    = panel_h * gap_frac;

    ax_r = axes('Position', [x0, y0_panel + fr_h + gap_h, panel_w, raster_h]);
    ax_f = axes('Position', [x0, y0_panel, panel_w, fr_h]);

    plot_raster_panel(ax_r, cd_cell, cfg, iC);
    plot_fr_panel(ax_f, cd_cell, cfg, iC);
end

% Panel 8: AHV heatmap
col = 4; row = 2;
x0 = margin_left + (col-1)*(panel_w + hspace);
y0_panel = margin_bottom + (nRows - row)*(panel_h + vspace);
ax_heat = axes('Position', [x0, y0_panel, panel_w, panel_h]);
plot_ahv_heatmap(ax_heat, cell_data, cfg);

end

%% ========================================================================
%  LOCAL FUNCTION: Plot raster panel
%  ========================================================================
function plot_raster_panel(ax, cd_cell, cfg, iC)
axes(ax); hold on;
nTrials = cd_cell.nTrials;

patch([0 cfg.laser_dur cfg.laser_dur 0], ...
      [0 0 nTrials+1 nTrials+1], ...
      [0.7 0.85 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

if ~isempty(cd_cell.raster_spk_times)
    plot(cd_cell.raster_spk_times, cd_cell.raster_trial_num, '.k', ...
        'MarkerSize', cfg.rasterTickSize);
end

if nTrials > 80;     ytks = [20 60 100 140];
elseif nTrials > 40; ytks = [10 30 50];
else;                ytks = [10 20];
end
ytks = ytks(ytks <= nTrials);

set(ax, 'XLim', cfg.window, 'YLim', [0.5 nTrials+0.5], ...
    'TickDir', 'out', 'FontSize', cfg.FontSize - 4, 'FontName', 'Arial', ...
    'XTickLabel', [], 'XTick', [-1 0 1 2], ...
    'YTick', ytks, 'Box', 'off');
if iC == 1 || iC == 5
    ylabel('Trial', 'FontSize', cfg.FontSize - 4, 'FontName', 'Arial');
end
if cfg.showTitles
    title(cd_cell.cell_label, 'FontSize', cfg.FontSize - 6, ...
        'Interpreter', 'none', 'FontWeight', 'normal');
end
end

%% ========================================================================
%  LOCAL FUNCTION: Plot FR PETH panel
%  ========================================================================
function plot_fr_panel(ax, cd_cell, cfg, iC)
axes(ax); hold on;

fr_ylim_max = max(cd_cell.fr_smooth) * 1.3;
if fr_ylim_max == 0; fr_ylim_max = 10; end

% Blue laser shading
patch([0 cfg.laser_dur cfg.laser_dur 0], ...
      [0 0 fr_ylim_max fr_ylim_max], ...
      [0.7 0.85 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Shuffle mean +/- SEM (black line + grey shading)
has_shuffle = ~isempty(cfg.shuff_results) && length(cfg.shuff_results) >= iC ...
              && ~isempty(cfg.shuff_results{iC});
if has_shuffle
    sr = cfg.shuff_results{iC};
    shuff_upper = sr.shuff_mean + sr.shuff_sem;
    shuff_lower = sr.shuff_mean - sr.shuff_sem;
    fill([sr.bin_centers fliplr(sr.bin_centers)], ...
         [shuff_upper fliplr(shuff_lower)], ...
         [0.75 0.75 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    plot(sr.bin_centers, sr.shuff_mean, 'k', 'LineWidth', 1.5);
end

% Real FR trace (red)
plot(cd_cell.bin_centers, cd_cell.fr_smooth, 'r', 'LineWidth', 2.5);

% Significance bars from shuffle (dark green p<.001, light green p<.05)
if has_shuffle
    green_y = fr_ylim_max * 0.95;
    dark_green  = [0.0 0.5 0.0];
    light_green = [0.4 0.8 0.3];
    bar_lw = 6;

    % Light green first (p < .05, including .001 bins — draw as background)
    sig_any = sr.sig_001 | sr.sig_05;
    if any(sig_any)
        draw_sig_bar(sr.bin_centers, sig_any, cfg.dt, green_y, light_green, bar_lw);
    end

    % Dark green on top (p < .001)
    if any(sr.sig_001)
        draw_sig_bar(sr.bin_centers, sr.sig_001, cfg.dt, green_y, dark_green, bar_lw);
    end
else
    % Placeholder solid green bar when no shuffle data
    green_y = fr_ylim_max * 0.95;
    plot([0 cfg.laser_dur], [green_y green_y], '-', ...
        'Color', [0.2 0.7 0.2], 'LineWidth', 4);
end

set(ax, 'XLim', cfg.window, 'YLim', [0 fr_ylim_max], ...
    'TickDir', 'out', 'FontSize', cfg.FontSize - 4, 'FontName', 'Arial', ...
    'XTick', [-1 0 1 2], 'Box', 'off');

% "Laser" label inside blue shading on bottom-left panel only
if iC == 5
    text(cfg.laser_dur / 2, fr_ylim_max * 0.15, 'Laser', ...
        'HorizontalAlignment', 'center', 'FontSize', cfg.FontSize - 4, ...
        'Color', [0.3 0.5 0.8], 'FontWeight', 'bold', 'FontName', 'Arial');
    xlabel('Time from laser onset (s)', 'FontSize', cfg.FontSize - 4, 'FontName', 'Arial');
end
if iC == 1 || iC == 5
    ylabel('FR (Hz)', 'FontSize', cfg.FontSize - 4, 'FontName', 'Arial');
end

% Cell number in top-left corner
text(cfg.window(1) + 0.15, fr_ylim_max * 0.88, num2str(iC), ...
    'FontSize', cfg.FontSize - 2, 'FontWeight', 'bold', ...
    'FontName', 'Arial', 'HorizontalAlignment', 'left');
end

%% ========================================================================
%  LOCAL FUNCTION: Plot AHV heatmap (version A, slot 8)
%  ========================================================================
function plot_ahv_heatmap(ax, cell_data, cfg)
axes(ax); hold on;

nCells_plot = length(cell_data);
bin_centers = cfg.AHV_binEdges(1:end-1) + diff(cfg.AHV_binEdges(1:2))/2;
nBins = length(bin_centers);
mid_bin = ceil(nBins / 2);  % center index

heat_matrix = nan(nCells_plot, nBins);
for iC = 1:nCells_plot
    if isempty(cell_data{iC}); continue; end
    cd_cell = cell_data{iC};
    tc_vals = cd_cell.AHV_tc_smooth;
    good = cd_cell.AHV_good_bins;
    valid_vals = tc_vals(good);
    if isempty(valid_vals) || range(valid_vals) == 0; continue; end

    % Determine preferred side: compare mean FR on negative vs positive AHV
    neg_bins = good & bin_centers < 0;
    pos_bins = good & bin_centers > 0;
    mean_neg = nanmean(tc_vals(neg_bins));
    mean_pos = nanmean(tc_vals(pos_bins));

    % Flip so preferred (higher FR) side is on the left (negative AHV)
    if mean_pos > mean_neg
        tc_oriented = fliplr(tc_vals);
        good_oriented = fliplr(good);
    else
        tc_oriented = tc_vals;
        good_oriented = good;
    end

    % Normalize to [0, 1]
    valid_oriented = tc_oriented(good_oriented);
    tc_norm = (tc_oriented - min(valid_oriented)) / range(valid_oriented);
    tc_norm(~good_oriented) = NaN;
    
    % Trim isolated edge bins (good bin, NaN gap, then more good bins).
    % Finds the first and last contiguous run from the inside out;
    % NaN-out anything outside that range.
    good_idx = find(~isnan(tc_norm));
    good_idx = good_idx(:);  % force column
    if length(good_idx) > 2
        diffs = diff(good_idx);
        gaps = find(diffs > 1);
        if ~isempty(gaps)
            % Keep only the largest contiguous block
            block_starts = [good_idx(1); good_idx(gaps(:)+1)];
            block_ends   = [good_idx(gaps(:)); good_idx(end)];
            block_lens   = block_ends - block_starts + 1;
            [~, best]    = max(block_lens);
            keep_start   = block_starts(best);
            keep_end     = block_ends(best);
            tc_norm(1:keep_start-1) = NaN;
            tc_norm(keep_end+1:end) = NaN;
        end
    end
    
    heat_matrix(iC, :) = tc_norm;
end

% No sorting — rows match PETH panel order (cell 1 = row 1, etc.)
heat_matrix_sorted = heat_matrix;

% Force white background behind transparent NaN pixels
set(ax, 'Color', 'w');
set(get(ax, 'Parent'), 'Color', 'w');

h_img = imagesc(bin_centers, 1:nCells_plot, heat_matrix_sorted);
set(h_img, 'AlphaData', ~isnan(heat_matrix_sorted));
set(ax, 'YDir', 'reverse', 'XLim', [-200 200], 'YLim', [0.5 nCells_plot+0.5], ...
    'TickDir', 'out', 'FontSize', cfg.FontSize - 4, ...
    'XTick', [-200 0 200], 'YTick', 1:nCells_plot, 'Box', 'off');
colormap(ax, flipud(gray));  % dark = high FR (matches manuscript Fig 2B)
set(ax, 'CLim', [0 1]);  % lock colormap range to normalized [0,1]
xlabel('AHV (deg/s)', 'FontSize', cfg.FontSize - 4);
ylabel('NPH unit #', 'FontSize', cfg.FontSize - 4);
end

%% ========================================================================
%  LOCAL FUNCTION: Draw significance bar as contiguous segments
%  ========================================================================
function draw_sig_bar(bin_centers, sig_mask, dt, y_val, color, lw)
% Draws contiguous runs of significant bins as single line segments
% to avoid gaps between adjacent bins.
sig_idx = find(sig_mask);
if isempty(sig_idx); return; end

% Find contiguous runs
breaks = find(diff(sig_idx) > 1);
run_starts = [sig_idx(1), sig_idx(breaks + 1)];
run_ends   = [sig_idx(breaks), sig_idx(end)];

for iR = 1:length(run_starts)
    x_start = bin_centers(run_starts(iR)) - dt/2;
    x_end   = bin_centers(run_ends(iR))   + dt/2;
    plot([x_start x_end], [y_val y_val], '-', 'Color', color, 'LineWidth', lw);
end
end

%% Version History
% 2026-03-18  v2.4
%                   - Real significance bars from compute_optotag_shuffle:
%                     dark green (p<.001), light green (p<.05)
%                   - Shuffle mean +/- SEM overlay (black line + grey shading)
%                   - Falls back to placeholder green bar if no shuff_results
%                   - New cfg field: .shuff_results (cell array)
% 2026-03-18  v2.3
%                   - Added spike_times and session_duration to cell_data
%                     (needed for compute_optotag_shuffle circular shuffle)
%                   - LoadSpikesQuick now returns spk_params for starttime
% 2026-03-18  v2.2
%                   - Three figure versions: A (PETHs only), B (TC insets +
%                     heatmap), C (PETHs + heatmap)
%                   - Heatmap: flipud(gray) so dark=high FR (matches Fig 2B)
%                   - Heatmap: YDir reverse (row 1 on top)
%                   - Heatmap: removed peak-bin sorting (rows match PETH order)
%                   - Heatmap: trim isolated edge bins (fixes cell 5 gaps)
%                   - Cell number label in top-left of each FR panel
%                   - ylabels only on left-column panels (iC 1 and 5)
% 2026-03-18  v2.1
%                   - Fixed heatmap NaN background (forced white axes + CLim)
%                   - Fixed "?" in Version B TC y-axis (FontName Arial on all axes)
% 2026-03-18  v2.0
%                   - Separated data computation from plotting
%                   - Added AHV TC computation per cell
%                   - Two figure versions: A (heatmap slot 8), B (TC inset)
%                   - Always saves .mat + .png per run
% 2026-03-18  v1.1
%                   - Titles optional, save fix, fewer ticks
% 2026-03-18  v1.0 Initial version.
