function [fig_handle, results, params] = plot_stationary_saccade_PETHs_ver_1_1(varargin)
% plot_stationary_saccade_PETHs_ver_1_1  Diagnostic: saccade PETHs during
% platform-stationary periods vs. all saccades.
%
%   For a single neuron, generates a 2-row x 2-column figure:
%     Row 1: All saccades (CCW / CW) — raster + FR overlay
%     Row 2: Stationary-only saccades (CCW / CW) — raster + FR overlay
%
%   This lets you compare the saccade-evoked burst pattern when the
%   platform is moving (VOR context) vs. stationary (spontaneous saccades).
%
%   INPUTS:
%       varargin    - name-value pairs (see defaults below)
%
%   KEY OPTIONS:
%       neuron_name         - cell name, e.g. 'M079-2020-10-28-TT07_2'
%       tfilelist_path      - path to NPH_all_tfilelist.mat
%       data_root           - root data path
%       save_root           - analysis output path
%       peth_window         - [pre post] in seconds (default [-0.2 0.2])
%       min_saccades        - minimum stationary saccades to attempt PETH
%       doSavePNG, doSaveFIG, doShow
%
%   OUTPUTS:
%       fig_handle  - figure handle (empty if insufficient stationary saccades)
%       results     - struct with all computed data
%       params      - struct of parameter values used
%
%   USAGE:
%       [fh, res, par] = plot_stationary_saccade_PETHs_ver_1_1('neuron_name', 'M079-2020-10-28-TT07_2', 'doShow', 1);
%
%   DEPENDENCIES:
%       LoadSpikesQuick, SelectTS, restrict, pushdir, popdir,
%       HD_GetSSN, process_varargin, tsd, nearest_idx3, MakeQfromS,
%       load_chenhao_saccades, save_fig_trio, find_cell_index

%% ========================================================================
%  DEFAULTS
%  ========================================================================
neuron_name         = 'M079-2020-10-28-TT07_2';

tfilelist_path      = 'C:\Jeff\U01\datatouse\NPH all tfilelist.mat';
chenhao_sacc_path   = 'C:\Jeff\U01\datatouse\saccade_ONSET_all_unit.mat';
data_root           = 'C:\Jeff\U01\datatouse';
save_root           = 'C:\Jeff\U01\datatouse\___analysis_output';

% Saccade PETH parameters
peth_window         = [-0.2 0.2];       % s
peth_dt             = 0.005;            % s, bin width for FR PETH
peth_smooth_win     = 5;                % bins, smoothdata gaussian window
peth_excessBounds   = 1;                % s, extra restrict margin

% Stationary saccade filter
min_saccades        = 3;                % minimum per direction to attempt PETH

% Figure parameters
FontSize            = 22;
TickFontSize        = 16;
TitleFontSize       = 24;
MarkerSize_raster   = 3;
LineWidth_FR        = 4;
color_all_CCW       = [0.30 0.69 0.29];     % green
color_all_CW        = [0.90 0.45 0.15];     % orange
color_stat_CCW      = [0.20 0.55 0.20];     % darker green
color_stat_CW       = [0.75 0.35 0.10];     % darker orange
rasterColor         = [0.25 0.25 0.25];
figWidth            = 1200;
figHeight           = 900;

% Output control
doSavePNG   = 1;
doSaveFIG   = 1;
doSaveEPS   = 1;
doShow      = 0;
doSave      = 1;
verbose     = 1;

%% Process varargin
process_varargin(varargin);

%% ========================================================================
%  PARAMS (analysis-relevant free parameters only)
%  ========================================================================
params = struct();
params.function_name        = 'plot_stationary_saccade_PETHs';
params.version              = '1.1';
params.date_run             = datestr(now, 'yyyy-mm-dd');
params.neuron_name          = neuron_name;
params.saccade_source       = 'Chenhao onset times (saccade_ONSET_all_unit.mat)';
params.peth_window          = peth_window;
params.peth_dt              = peth_dt;
params.peth_smooth_win      = peth_smooth_win;
params.peth_excessBounds    = peth_excessBounds;
params.min_saccades         = min_saccades;

%% ========================================================================
%  FIND SESSION PATH FROM TFILELIST
%  ========================================================================
tfl = load(tfilelist_path, 'NPH_all_tfilelist');
NPH_all_tfilelist = tfl.NPH_all_tfilelist;

sess_path = '';
for kk = 1:length(NPH_all_tfilelist)
    [kk_folder, kk_name] = fileparts(NPH_all_tfilelist{kk});
    if strcmp(kk_name, neuron_name)
        sess_path = kk_folder;
        break
    end
end
if isempty(sess_path)
    error('Could not find neuron %s in tfilelist.', neuron_name);
end

if verbose; fprintf('\n=== Stationary Saccade PETH Diagnostic ===\n'); end
if verbose; fprintf('Neuron: %s\n', neuron_name); end
if verbose; fprintf('Session: %s\n', sess_path); end

%% ========================================================================
%  LOAD DATA
%  ========================================================================
pushdir(sess_path);
SSN = HD_GetSSN;

% --- Spikes ---
[S_all, spike_cell_names, spk_params] = LoadSpikesQuick();

cell_idx = find_cell_index(spike_cell_names, neuron_name, SSN);
if isempty(cell_idx)
    popdir;
    error('Could not match %s to loaded spike trains. Names: %s', ...
        neuron_name, strjoin(spike_cell_names, ', '));
end
myCell = SelectTS([], S_all, cell_idx);

% --- Saccade onset times (Chenhao) ---
[temporal_all, nasal_all] = load_chenhao_saccades(SSN, ...
    'chenhao_sacc_path', chenhao_sacc_path);
if verbose
    fprintf('  Chenhao saccade onsets: %d CCW, %d CW\n', ...
        length(temporal_all), length(nasal_all));
end

% --- Stationary times ---
stat_file = strcat(SSN, '-AHV_StationaryTimes.mat');
if ~exist(stat_file, 'file')
    popdir;
    error('Stationary times file not found: %s', stat_file);
end
stat_data = load(stat_file, 'STtstart', 'STtend');
STtstart = stat_data.STtstart;
STtend   = stat_data.STtend;

%% ========================================================================
%  ISOLATE STATIONARY SACCADES
%  (same logic as isolateStationarySaccades.m)
%  ========================================================================
temporal_stat_idx = false(size(temporal_all));
nasal_stat_idx    = false(size(nasal_all));

for iST = 1:length(STtstart)
    temporal_stat_idx = temporal_stat_idx | ...
        (temporal_all > STtstart(iST) & temporal_all < STtend(iST));
    nasal_stat_idx = nasal_stat_idx | ...
        (nasal_all > STtstart(iST) & nasal_all < STtend(iST));
end

temporal_stat = temporal_all(temporal_stat_idx);
nasal_stat    = nasal_all(nasal_stat_idx);

% Also get moving-only saccades for reference
temporal_moving = temporal_all(~temporal_stat_idx);
nasal_moving    = nasal_all(~nasal_stat_idx);

%% Report counts
nTemporal_all  = length(temporal_all);
nNasal_all     = length(nasal_all);
nTemporal_stat = length(temporal_stat);
nNasal_stat    = length(nasal_stat);
nTemporal_move = length(temporal_moving);
nNasal_move    = length(nasal_moving);

if verbose
    fprintf('\n  Saccade counts:\n');
    fprintf('    CCW (temporal):  %d total  |  %d moving  |  %d stationary\n', ...
        nTemporal_all, nTemporal_move, nTemporal_stat);
    fprintf('    CW  (nasal):     %d total  |  %d moving  |  %d stationary\n', ...
        nNasal_all, nNasal_move, nNasal_stat);
    fprintf('    Stationary periods: %d\n', length(STtstart));
    
    % Total stationary duration
    stat_dur_total = sum(STtend - STtstart);
    fprintf('    Total stationary duration: %.1f s\n', stat_dur_total);
end

% Store in results
results = struct();
results.neuron_name         = neuron_name;
results.SSN                 = SSN;
results.sess_path           = sess_path;
results.nTemporal_all       = nTemporal_all;
results.nNasal_all          = nNasal_all;
results.nTemporal_stat      = nTemporal_stat;
results.nNasal_stat         = nNasal_stat;
results.nTemporal_move      = nTemporal_move;
results.nNasal_move         = nNasal_move;
results.n_stationary_periods = length(STtstart);
results.stationary_dur_total = sum(STtend - STtstart);
results.temporal_stat_times  = temporal_stat;
results.nasal_stat_times     = nasal_stat;
results.temporal_all_times   = temporal_all;
results.nasal_all_times      = nasal_all;

%% ========================================================================
%  CHECK MINIMUM SACCADE COUNTS
%  ========================================================================
has_enough_stat = (nTemporal_stat >= min_saccades) || (nNasal_stat >= min_saccades);

if ~has_enough_stat
    fprintf('\n  *** Insufficient stationary saccades (need >= %d in at least one direction). ***\n', min_saccades);
    fprintf('  *** Skipping figure. ***\n');
    fig_handle = [];
    results.has_enough_stationary = false;
    popdir;
    return
end
results.has_enough_stationary = true;

%% ========================================================================
%  COMPUTE RASTERS AND PETHs
%  ========================================================================
peth_bins    = linspace(peth_window(1), peth_window(2), diff(peth_window)/peth_dt + 1);
peth_centers = peth_bins(1:end-1) + diff(peth_bins) / 2;

% --- All saccades ---
[spk_rel_T_all, evt_idx_T_all] = compute_raster_data(myCell, temporal_all, peth_window, peth_excessBounds);
[m_T_all, ~] = histcounts(spk_rel_T_all, peth_bins);
fr_T_all = m_T_all / peth_dt / max(nTemporal_all, 1);
fr_T_all_smooth = smoothdata(movmean(fr_T_all, peth_smooth_win), 'gaussian', peth_smooth_win);

[spk_rel_N_all, evt_idx_N_all] = compute_raster_data(myCell, nasal_all, peth_window, peth_excessBounds);
[m_N_all, ~] = histcounts(spk_rel_N_all, peth_bins);
fr_N_all = m_N_all / peth_dt / max(nNasal_all, 1);
fr_N_all_smooth = smoothdata(movmean(fr_N_all, peth_smooth_win), 'gaussian', peth_smooth_win);

% --- Stationary saccades ---
if nTemporal_stat >= min_saccades
    [spk_rel_T_stat, evt_idx_T_stat] = compute_raster_data(myCell, temporal_stat, peth_window, peth_excessBounds);
    [m_T_stat, ~] = histcounts(spk_rel_T_stat, peth_bins);
    fr_T_stat = m_T_stat / peth_dt / max(nTemporal_stat, 1);
    fr_T_stat_smooth = smoothdata(movmean(fr_T_stat, peth_smooth_win), 'gaussian', peth_smooth_win);
    has_T_stat = true;
else
    spk_rel_T_stat = []; evt_idx_T_stat = [];
    fr_T_stat_smooth = nan(size(peth_centers));
    has_T_stat = false;
    if verbose; fprintf('  Not enough stationary CCW saccades (%d < %d)\n', nTemporal_stat, min_saccades); end
end

if nNasal_stat >= min_saccades
    [spk_rel_N_stat, evt_idx_N_stat] = compute_raster_data(myCell, nasal_stat, peth_window, peth_excessBounds);
    [m_N_stat, ~] = histcounts(spk_rel_N_stat, peth_bins);
    fr_N_stat = m_N_stat / peth_dt / max(nNasal_stat, 1);
    fr_N_stat_smooth = smoothdata(movmean(fr_N_stat, peth_smooth_win), 'gaussian', peth_smooth_win);
    has_N_stat = true;
else
    spk_rel_N_stat = []; evt_idx_N_stat = [];
    fr_N_stat_smooth = nan(size(peth_centers));
    has_N_stat = false;
    if verbose; fprintf('  Not enough stationary CW saccades (%d < %d)\n', nNasal_stat, min_saccades); end
end

% Store PETH data
results.peth_centers        = peth_centers;
results.fr_T_all_smooth     = fr_T_all_smooth;
results.fr_N_all_smooth     = fr_N_all_smooth;
results.fr_T_stat_smooth    = fr_T_stat_smooth;
results.fr_N_stat_smooth    = fr_N_stat_smooth;
results.has_T_stat          = has_T_stat;
results.has_N_stat          = has_N_stat;

popdir;

%% ========================================================================
%  CREATE FIGURE
%  ========================================================================
if doShow
    fig_handle = figure('Position', [50 50 figWidth figHeight]);
else
    fig_handle = figure('Position', [50 50 figWidth figHeight], 'Visible', 'off');
end

% Store metadata
fig_handle.UserData.created   = datestr(now, 'yyyy-mm-dd HH:MM:SS');
fig_handle.UserData.function  = 'plot_stationary_saccade_PETHs_ver_1_1';
fig_handle.UserData.neuron    = neuron_name;

% Layout: 2 rows x 2 columns
% Row 1: All saccades (CCW, CW)
% Row 2: Stationary saccades (CCW, CW)

row_labels = {'All Saccades', 'Stationary Only'};
col_labels = {'CCW (Temporal)', 'CW (Nasal)'};

% Raster + FR data for each panel {row, col}
panel_data = {
    % Row 1: All
    struct('spk_rel', spk_rel_T_all, 'evt_idx', evt_idx_T_all, ...
           'fr_smooth', fr_T_all_smooth, 'nEvents', nTemporal_all, ...
           'fr_color', color_all_CCW, 'has_data', true), ...
    struct('spk_rel', spk_rel_N_all, 'evt_idx', evt_idx_N_all, ...
           'fr_smooth', fr_N_all_smooth, 'nEvents', nNasal_all, ...
           'fr_color', color_all_CW, 'has_data', true);
    % Row 2: Stationary
    struct('spk_rel', spk_rel_T_stat, 'evt_idx', evt_idx_T_stat, ...
           'fr_smooth', fr_T_stat_smooth, 'nEvents', nTemporal_stat, ...
           'fr_color', color_stat_CCW, 'has_data', has_T_stat), ...
    struct('spk_rel', spk_rel_N_stat, 'evt_idx', evt_idx_N_stat, ...
           'fr_smooth', fr_N_stat_smooth, 'nEvents', nNasal_stat, ...
           'fr_color', color_stat_CW, 'has_data', has_N_stat);
};

% Match FR y-axes across all 4 panels for direct comparison
all_fr_max = 0;
for iR = 1:2
    for iC = 1:2
        pd = panel_data{iR, iC};
        if pd.has_data && ~all(isnan(pd.fr_smooth))
            all_fr_max = max(all_fr_max, max(pd.fr_smooth));
        end
    end
end
fr_ylim_upper = ceil(all_fr_max / 10) * 10;  % round up to nearest 10
if fr_ylim_upper == 0; fr_ylim_upper = 50; end

% Axes positioning
left_margin     = 0.10;
right_margin    = 0.08;
top_margin      = 0.12;
bottom_margin   = 0.10;
row_gap         = 0.08;
col_gap         = 0.06;

total_height    = 1 - top_margin - bottom_margin;
row_height      = (total_height - row_gap) / 2;
total_width     = 1 - left_margin - right_margin;
col_width       = (total_width - col_gap) / 2;

for iR = 1:2
    for iC = 1:2
        x_pos = left_margin + (iC - 1) * (col_width + col_gap);
        y_pos = 1 - top_margin - iR * row_height - (iR - 1) * row_gap;

        ax = axes('Position', [x_pos, y_pos, col_width, row_height]);
        hold(ax, 'on');

        pd = panel_data{iR, iC};

        if pd.has_data && ~isempty(pd.spk_rel)
            % --- Raster on left y-axis ---
            yyaxis(ax, 'left');
            plot(ax, pd.spk_rel, pd.evt_idx + 0.5, '.', ...
                'Color', rasterColor, 'MarkerSize', MarkerSize_raster);
            ylim(ax, [0 pd.nEvents + 1]);
            ax.YAxis(1).Color = [0 0 0];

            if iC == 1
                ylabel(ax, 'Trial', 'FontSize', FontSize);
            else
                set(ax.YAxis(1), 'TickLabels', []);
            end

            % --- FR on right y-axis ---
            yyaxis(ax, 'right');
            plot(ax, peth_centers, pd.fr_smooth, '-', ...
                'Color', pd.fr_color, 'LineWidth', LineWidth_FR);
            ylim(ax, [0 fr_ylim_upper]);
            ax.YAxis(2).Color = pd.fr_color;

            if iC == 2
                ylabel(ax, 'FR (Hz)', 'FontSize', FontSize);
            else
                set(ax.YAxis(2), 'TickLabels', []);
            end

        elseif pd.has_data == false
            % Not enough saccades
            text(ax, 0.5, 0.5, sprintf('n = %d\n(< %d)', pd.nEvents, min_saccades), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'FontSize', FontSize, 'Color', [0.6 0.6 0.6]);
        end

        % Zero line
        xline(ax, 0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

        % Axis formatting
        xlim(ax, peth_window);
        ax.XAxis.FontSize = TickFontSize;
        if iR == 2
            xlabel(ax, 'Time (s)', 'FontSize', FontSize);
        else
            set(ax, 'XTickLabel', []);
        end

        % Panel title
        if iR == 1
            title_str = sprintf('%s  (n = %d)', col_labels{iC}, pd.nEvents);
        else
            title_str = sprintf('%s  (n = %d)', col_labels{iC}, pd.nEvents);
        end
        title(ax, title_str, 'FontSize', TitleFontSize, 'FontWeight', 'bold');

        box(ax, 'off');
    end
end

% Row labels (left side, rotated via invisible axes + text)
ax_label = axes('Position', [0 0 1 1], 'Visible', 'off');
for iR = 1:2
    y_center = 1 - top_margin - (iR - 0.5) * row_height - (iR - 1) * row_gap;
    text(ax_label, 0.02, y_center, row_labels{iR}, ...
        'FontSize', FontSize, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Rotation', 90, 'Units', 'normalized');
end

% Supertitle
sgtitle_str = sprintf('%s — Stationary vs All Saccades', strrep(neuron_name, '_', '\_'));
sgtitle(sgtitle_str, 'FontSize', TitleFontSize + 2, 'FontWeight', 'bold');

%% ========================================================================
%  SAVE
%  ========================================================================
if doSave
    % Create output folder
    date_str = datestr(now, 'yyyy-mm-dd');
    base_folder_name = sprintf('stationary_saccade_PETHs_%s', date_str);
    run_num = 1;
    out_folder = fullfile(save_root, sprintf('%s_run%d', base_folder_name, run_num));
    while exist(out_folder, 'dir')
        run_num = run_num + 1;
        out_folder = fullfile(save_root, sprintf('%s_run%d', base_folder_name, run_num));
    end
    mkdir(out_folder);

    safe_name = strrep(neuron_name, '-', '_');
    filepath_base = fullfile(out_folder, sprintf('stationary_PETH_%s', safe_name));
    save_fig_trio(fig_handle, filepath_base, ...
        'doSavePNG', doSavePNG, 'doSaveFIG', doSaveFIG, 'doSaveEPS', doSaveEPS, 'doShow', doShow);

    % Save results and params
    save(fullfile(out_folder, sprintf('stationary_PETH_results_%s.mat', safe_name)), ...
        'results', 'params', '-v7');
    if verbose; fprintf('  Saved to: %s\n', out_folder); end

    results.output_folder = out_folder;
end

if verbose; fprintf('\n  Done.\n'); end

end

%% ========================================================================
%  LOCAL HELPER: compute raster data
%  (same logic as plot_fig4B_ver_2_2)
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

%% Version History
% 2026-04-05  v1.1  Refactor:
%                    - Switched to Chenhao's saccade onset times via
%                      load_chenhao_saccades.m (Jeff's saccades-edited.mat
%                      times are now deprecated for saccade timestamps)
%                    - Cell matching via find_cell_index.m (shared utility)
%                    - Figure saving via save_fig_trio.m (shared utility,
%                      includes savefig visibility bug workaround)
%                    - doSaveEPS default changed to 1 (was 0)
%                    - .mat save now uses -v7 flag
%                    - Removed unused has_eye_traces variable
% 2026-04-03  v1.0  Initial version. Single-neuron diagnostic comparing
%                    saccade PETHs during platform-stationary vs. all
%                    periods. Raster + FR overlay, 2x2 layout (all/stat x
%                    CCW/CW). Uses isolateStationarySaccades logic inline.
