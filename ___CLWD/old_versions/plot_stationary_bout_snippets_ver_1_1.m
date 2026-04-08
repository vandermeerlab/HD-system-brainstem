function [fig_handle, results, params] = plot_stationary_bout_snippets_ver_1_1(varargin)
% plot_stationary_bout_snippets_ver_1_1  Visualize stationary saccade bouts
% with eye velocity, firing rate, and saccade onset markers.
%
%   For a single neuron, finds the best platform-stationary periods
%   (ranked by total saccade count), and plots time-series snippets
%   showing eye velocity + smoothed firing rate + saccade onset markers
%   + spike raster ticks.
%
%   Uses Chenhao's revised saccade onset times (saccade_ONSET_all_unit.mat).
%
%   INPUTS (name-value pairs via varargin):
%       neuron_name         - cell name (e.g., 'M402-2022-11-02-2-TT02_3')
%       neuron_index        - index into NPH_all_tfilelist (alternative)
%       tfilelist_path      - path to NPH all tfilelist.mat
%       chenhao_sacc_path   - path to saccade_ONSET_all_unit.mat
%       n_bouts             - how many best bouts to show (default: 3)
%       snippet_duration    - max seconds per panel (default: 15)
%       min_bout_saccades   - minimum saccades per bout to qualify (default: 3)
%       fr_bin_dt           - firing rate bin width in s (default: 0.01)
%       fr_smooth_win       - Gaussian smoothing window in bins (default: 15)
%       ev_padding          - fraction of range to pad eye vel y-axis (default: 0.15)
%       doSavePNG, doSaveFIG, doSaveEPS, doShow
%
%   OUTPUTS:
%       fig_handle  - figure handle ([] if insufficient data)
%       results     - struct with bout info, saccade counts, etc.
%       params      - struct of parameter values used
%
%   USAGE:
%       [fh, res, par] = plot_stationary_bout_snippets_ver_1_1(...
%           'neuron_name', 'M402-2022-11-02-2-TT02_3', 'doShow', 1);
%
%   DEPENDENCIES:
%       LoadSpikesQuick, SelectTS, restrict, pushdir, popdir,
%       HD_GetSSN, process_varargin

%% ========================================================================
%  DEFAULTS
%  ========================================================================
neuron_name         = '';
neuron_index        = [];

tfilelist_path      = 'C:\Jeff\U01\datatouse\NPH all tfilelist.mat';
chenhao_sacc_path   = 'C:\Jeff\U01\datatouse\saccade_ONSET_all_unit.mat';
save_root           = 'C:\Jeff\U01\datatouse\___analysis_output';

% Bout selection
n_bouts             = 3;
snippet_duration    = 15;       % max seconds per snippet window
min_bout_saccades   = 3;        % minimum saccades in a bout to qualify

% Firing rate computation
fr_bin_dt           = 0.01;     % s
fr_smooth_win       = 15;       % bins for Gaussian smoothing

% Eye velocity display
ev_padding          = 0.15;     % fraction of data range to pad y-axis

% Figure parameters
FontSize            = 20;
TickFontSize        = 16;
TitleFontSize       = 18;
LineWidth_EV        = 1.2;
LineWidth_FR        = 3;
LineWidth_sacc      = 2.0;
alpha_sacc          = 0.55;
color_CCW_sacc      = [0.30 0.69 0.29];   % green — temporal / CCW
color_CW_sacc       = [0.90 0.45 0.15];   % orange — nasal / CW
color_FR            = [0.15 0.35 0.70];    % blue — firing rate
color_EV            = [0.35 0.35 0.35];    % grey — eye velocity
color_spk_tick      = [0.10 0.10 0.10];    % near-black — spike ticks
figWidth            = 1500;

% Output control
doSavePNG   = 1;
doSaveFIG   = 1;
doSaveEPS   = 1;
doShow      = 0;
verbose     = 1;

%% Process varargin
process_varargin(varargin);

%% ========================================================================
%  PARAMS
%  ========================================================================
params = struct();
params.function_name        = 'plot_stationary_bout_snippets';
params.version              = '1.1';
params.date_run             = datestr(now, 'yyyy-mm-dd HH:MM:SS');
params.n_bouts              = n_bouts;
params.snippet_duration     = snippet_duration;
params.min_bout_saccades    = min_bout_saccades;
params.fr_bin_dt            = fr_bin_dt;
params.fr_smooth_win        = fr_smooth_win;
params.ev_padding           = ev_padding;
params.saccade_source       = 'Chenhao onset times (saccade_ONSET_all_unit.mat)';

%% ========================================================================
%  RESOLVE NEURON
%  ========================================================================
tfl = load(tfilelist_path, 'NPH_all_tfilelist');
NPH_all_tfilelist = tfl.NPH_all_tfilelist(:);

if ~isempty(neuron_index) && isempty(neuron_name)
    [~, neuron_name] = fileparts(NPH_all_tfilelist{neuron_index});
end

if isempty(neuron_name)
    error('Must specify neuron_name or neuron_index.');
end

% Find session path
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

params.neuron_name = neuron_name;
if verbose; fprintf('\n=== Stationary Bout Snippets (v1.1) ===\n'); end
if verbose; fprintf('Neuron: %s\n', neuron_name); end

%% ========================================================================
%  LOAD CHENHAO SACCADE ONSET TIMES
%  ========================================================================
chenhao_data = load(chenhao_sacc_path, 'NPH_list', 'nasalSaccades_ONSET', 'temporalSaccades_ONSET');
chenhao_sessions = cell(1, length(chenhao_data.NPH_list));
for ii = 1:length(chenhao_data.NPH_list)
    chenhao_sessions{ii} = chenhao_data.NPH_list{ii};
end

%% ========================================================================
%  LOAD SESSION DATA
%  ========================================================================
pushdir(sess_path);
SSN = HD_GetSSN;
params.SSN = SSN;
if verbose; fprintf('Session: %s\n', SSN); end

% --- Find Chenhao session index ---
chenhao_idx = [];
for ii = 1:length(chenhao_sessions)
    if strcmp(chenhao_sessions{ii}, SSN)
        chenhao_idx = ii;
        break
    end
end
if isempty(chenhao_idx)
    popdir;
    error('Session %s not found in Chenhao saccade onset file.', SSN);
end

temporal_all = chenhao_data.temporalSaccades_ONSET{chenhao_idx}(:);  % CCW
nasal_all    = chenhao_data.nasalSaccades_ONSET{chenhao_idx}(:);     % CW
temporal_all = temporal_all(~isnan(temporal_all));
nasal_all    = nasal_all(~isnan(nasal_all));

if verbose
    fprintf('  Chenhao saccade onsets: %d CCW (temporal), %d CW (nasal)\n', ...
        length(temporal_all), length(nasal_all));
end

% --- Load spikes ---
[S_all, spike_cell_names, ~] = LoadSpikesQuick();
cell_idx = find_cell_index(spike_cell_names, neuron_name, SSN);
if isempty(cell_idx)
    popdir;
    error('Could not match %s to loaded spike trains.', neuron_name);
end
myCell = SelectTS([], S_all, cell_idx);

% --- Load eye velocity ---
sacc_file = strcat(SSN, '-saccades-edited.mat');
has_eye_data = false;
eye_vel_tsd = [];
if exist(sacc_file, 'file')
    try
        eye_fields = load(sacc_file, 'diffHdeg', 'diffH');
        if isfield(eye_fields, 'diffHdeg')
            eye_vel_tsd = eye_fields.diffHdeg;
            has_eye_data = true;
        elseif isfield(eye_fields, 'diffH')
            eye_vel_tsd = eye_fields.diffH;
            has_eye_data = true;
        end
    catch ME_eye
        if verbose; fprintf('  Could not load eye data: %s\n', ME_eye.message); end
    end
end

% --- Load stationary times ---
stat_file = strcat(SSN, '-AHV_StationaryTimes.mat');
if ~exist(stat_file, 'file')
    popdir;
    error('Stationary times file not found for %s.', SSN);
end
stat_data = load(stat_file, 'STtstart', 'STtend');
STtstart = stat_data.STtstart;
STtend   = stat_data.STtend;

% --- Identify stationary saccades ---
t_stat_mask = false(size(temporal_all));
n_stat_mask = false(size(nasal_all));
for iST = 1:length(STtstart)
    t_stat_mask = t_stat_mask | (temporal_all > STtstart(iST) & temporal_all < STtend(iST));
    n_stat_mask = n_stat_mask | (nasal_all > STtstart(iST) & nasal_all < STtend(iST));
end
temporal_stat = temporal_all(t_stat_mask);
nasal_stat    = nasal_all(n_stat_mask);

if verbose
    fprintf('  Stationary saccades: %d CCW, %d CW (%.1f s total stationary)\n', ...
        length(temporal_stat), length(nasal_stat), sum(STtend - STtstart));
end

popdir;

%% ========================================================================
%  RANK STATIONARY BOUTS BY SACCADE COUNT
%  ========================================================================
n_periods = length(STtstart);
sacc_per_period = zeros(n_periods, 1);

for iP = 1:n_periods
    nT_p = sum(temporal_stat > STtstart(iP) & temporal_stat < STtend(iP));
    nN_p = sum(nasal_stat > STtstart(iP) & nasal_stat < STtend(iP));
    sacc_per_period(iP) = nT_p + nN_p;
end

[~, sort_idx] = sort(sacc_per_period, 'descend');

% Pick top bouts with enough saccades
bout_indices = [];
for ii = 1:length(sort_idx)
    if sacc_per_period(sort_idx(ii)) >= min_bout_saccades
        bout_indices(end+1) = sort_idx(ii); %#ok<AGROW>
    end
    if length(bout_indices) >= n_bouts; break; end
end

if isempty(bout_indices)
    if verbose; fprintf('  No stationary bouts with >= %d saccades. Skipping.\n', min_bout_saccades); end
    fig_handle = []; results = struct(); return
end
n_actual_bouts = length(bout_indices);

%% ========================================================================
%  STORE RESULTS
%  ========================================================================
results = struct();
results.neuron_name     = neuron_name;
results.SSN             = SSN;
results.n_CCW_stat      = length(temporal_stat);
results.n_CW_stat       = length(nasal_stat);
results.n_stat_periods  = n_periods;
results.stat_dur_total  = sum(STtend - STtstart);
results.bout_indices    = bout_indices;
results.sacc_per_bout   = sacc_per_period(bout_indices);
results.saccade_source  = 'Chenhao onset times';

%% ========================================================================
%  PRECOMPUTE DATA FOR EACH BOUT (for consistent y-axis scaling)
%  ========================================================================
bout_data = cell(n_actual_bouts, 1);

if has_eye_data
    ev_tvec = eye_vel_tsd.tvec(:);
    ev_data = eye_vel_tsd.data(:);
end

for iB = 1:n_actual_bouts
    bd = struct();
    pIdx = bout_indices(iB);
    bd.t_start_full = STtstart(pIdx);
    bd.t_end_full   = STtend(pIdx);
    bd.period_dur   = bd.t_end_full - bd.t_start_full;

    % Get saccades in this period
    bd.sacc_T = temporal_stat(temporal_stat > bd.t_start_full & temporal_stat < bd.t_end_full);
    bd.sacc_N = nasal_stat(nasal_stat > bd.t_start_full & nasal_stat < bd.t_end_full);

    % Determine display window
    t_win_start = bd.t_start_full;
    t_win_end   = bd.t_end_full;
    if bd.period_dur > snippet_duration
        all_sacc_in = sort([bd.sacc_T; bd.sacc_N]);
        if ~isempty(all_sacc_in)
            mid_sacc = median(all_sacc_in);
            t_win_start = max(bd.t_start_full, mid_sacc - snippet_duration/2);
            t_win_end   = min(bd.t_end_full, t_win_start + snippet_duration);
        end
    end
    bd.t_start = t_win_start;
    bd.t_end   = t_win_end;

    % Re-filter saccades to display window
    bd.sacc_T_win = bd.sacc_T(bd.sacc_T > t_win_start & bd.sacc_T < t_win_end);
    bd.sacc_N_win = bd.sacc_N(bd.sacc_N > t_win_start & bd.sacc_N < t_win_end);

    % Eye velocity in window
    bd.has_ev = false;
    if has_eye_data
        ev_mask = ev_tvec >= t_win_start & ev_tvec <= t_win_end;
        if any(ev_mask)
            bd.ev_t = ev_tvec(ev_mask);
            bd.ev_d = ev_data(ev_mask);
            bd.has_ev = true;
        end
    end

    % Firing rate in window
    fr_bin_edges = t_win_start:fr_bin_dt:t_win_end;
    bd.has_fr = false;
    if length(fr_bin_edges) > 2
        spk_in_win = myCell.t{1}(myCell.t{1} >= t_win_start & myCell.t{1} <= t_win_end);
        fr_bin_counts = histcounts(spk_in_win, fr_bin_edges);
        bd.fr_raw = fr_bin_counts / fr_bin_dt;
        bd.fr_smooth = smoothdata(bd.fr_raw, 'gaussian', fr_smooth_win);
        bd.fr_smooth(bd.fr_smooth < 0) = 0;  % prevent negative FR from smoothing
        bd.fr_centers = fr_bin_edges(1:end-1) + fr_bin_dt/2;
        bd.spk_times = spk_in_win;
        bd.has_fr = true;
    end

    bout_data{iB} = bd;
end

% --- Compute shared FR y-axis across all bouts ---
all_fr_max = 0;
for iB = 1:n_actual_bouts
    bd = bout_data{iB};
    if bd.has_fr
        all_fr_max = max(all_fr_max, max(bd.fr_smooth));
    end
end
fr_ylim_upper = ceil(all_fr_max / 10) * 10;
if fr_ylim_upper < 10; fr_ylim_upper = 10; end

%% ========================================================================
%  FIGURE
%  ========================================================================
panel_height = 300;
fig_h = panel_height * n_actual_bouts + 180;
if doShow
    fig_handle = figure('Position', [50 50 figWidth fig_h]);
else
    fig_handle = figure('Position', [50 50 figWidth fig_h], 'Visible', 'off');
end
fig_handle.UserData.created        = datestr(now, 'yyyy-mm-dd HH:MM:SS');
fig_handle.UserData.function       = 'plot_stationary_bout_snippets_ver_1_1';
fig_handle.UserData.neuron         = neuron_name;
fig_handle.UserData.saccade_source = 'Chenhao onset times';

% Layout
lm = 0.08; rm = 0.08; tm = 0.07; bm = 0.07;
gap = 0.05;
usable_h = 1 - tm - bm;
ph = (usable_h - gap * (n_actual_bouts - 1)) / n_actual_bouts;

for iB = 1:n_actual_bouts
    bd = bout_data{iB};
    ypos = 1 - tm - iB * ph - (iB - 1) * gap;
    ax = axes('Position', [lm ypos (1 - lm - rm) ph]); hold(ax, 'on'); %#ok<LAXES>

    % --- Eye velocity (left y-axis) ---
    yyaxis(ax, 'left');
    if bd.has_ev
        plot(ax, bd.ev_t, bd.ev_d, '-', 'Color', color_EV, 'LineWidth', LineWidth_EV);
        % Auto-scale: use data range + padding
        ev_range = max(bd.ev_d) - min(bd.ev_d);
        if ev_range < 10; ev_range = 100; end
        ev_lo = min(bd.ev_d) - ev_padding * ev_range;
        ev_hi = max(bd.ev_d) + ev_padding * ev_range;
        ylim(ax, [ev_lo ev_hi]);
    end
    ylabel(ax, 'Eye Vel (\circ/s)', 'FontSize', FontSize);
    ax.YAxis(1).Color = [0 0 0];
    ax.YAxis(1).FontSize = TickFontSize;

    % --- Saccade onset markers (xline spans full panel automatically) ---
    for iS = 1:length(bd.sacc_T_win)
        xline(ax, bd.sacc_T_win(iS), '-', 'Color', [color_CCW_sacc alpha_sacc], ...
            'LineWidth', LineWidth_sacc);
    end
    for iS = 1:length(bd.sacc_N_win)
        xline(ax, bd.sacc_N_win(iS), '-', 'Color', [color_CW_sacc alpha_sacc], ...
            'LineWidth', LineWidth_sacc);
    end

    % --- Spike raster ticks along top of panel ---
    if bd.has_fr && ~isempty(bd.spk_times)
        yl = ylim(ax);
        tick_y_base = yl(2) - 0.03 * diff(yl);
        tick_y_top  = yl(2);
        for iSpk = 1:length(bd.spk_times)
            plot(ax, [bd.spk_times(iSpk) bd.spk_times(iSpk)], ...
                [tick_y_base tick_y_top], '-', 'Color', [color_spk_tick 0.35], 'LineWidth', 0.5);
        end
    end

    % --- Firing rate (right y-axis) ---
    yyaxis(ax, 'right');
    if bd.has_fr
        plot(ax, bd.fr_centers, bd.fr_smooth, '-', ...
            'Color', color_FR, 'LineWidth', LineWidth_FR);
        ylim(ax, [0 fr_ylim_upper]);
    end
    ylabel(ax, 'FR (Hz)', 'FontSize', FontSize);
    ax.YAxis(2).Color = color_FR;
    ax.YAxis(2).FontSize = TickFontSize;

    % --- Formatting ---
    xlim(ax, [bd.t_start bd.t_end]);
    ax.FontSize = TickFontSize;
    ax.Box = 'off';

    title_str = sprintf('Bout %d  —  %.1fs  |  %d CCW + %d CW saccades', ...
        iB, bd.period_dur, length(bd.sacc_T_win), length(bd.sacc_N_win));
    title(ax, title_str, 'FontSize', TitleFontSize, 'FontWeight', 'bold');

    if iB == n_actual_bouts
        xlabel(ax, 'Time (s)', 'FontSize', FontSize);
    else
        set(ax, 'XTickLabel', []);
    end
end

% --- Supertitle ---
sgtitle(sprintf('%s — Stationary Saccade Bouts (Chenhao onsets)', ...
    strrep(neuron_name, '_', '\_')), 'FontSize', TitleFontSize + 4, 'FontWeight', 'bold');

%% ========================================================================
%  SAVE
%  ========================================================================
date_str = datestr(now, 'yyyy-mm-dd');
run_num = 1;
out_folder = fullfile(save_root, sprintf('bout_snippets_%s_run%d', date_str, run_num));
while exist(out_folder, 'dir')
    run_num = run_num + 1;
    out_folder = fullfile(save_root, sprintf('bout_snippets_%s_run%d', date_str, run_num));
end
mkdir(out_folder);

safe_name = strrep(neuron_name, '-', '_');
filepath_base = fullfile(out_folder, sprintf('snippet_%s', safe_name));
if isvalid(fig_handle)
    if doSavePNG; print(fig_handle, [filepath_base '.png'], '-dpng', '-r600'); end
    if doSaveFIG
        try savefig(fig_handle, [filepath_base '.fig']);
        catch; set(fig_handle, 'Visible', 'on'); savefig(fig_handle, [filepath_base '.fig']);
            if ~doShow; set(fig_handle, 'Visible', 'off'); end
        end
    end
    if doSaveEPS; print(fig_handle, [filepath_base '.eps'], '-depsc', '-painters'); end
else
    warning('Figure handle invalid — skipping save.');
end

save(fullfile(out_folder, 'snippet_results.mat'), 'results', 'params', '-v7');

if verbose; fprintf('\n  Saved to: %s\n', out_folder); end
if ~doShow; close(fig_handle); end
if verbose; fprintf('  Done.\n'); end

end

%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================

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

%% Version History
% 2026-04-04  v1.1  Improvements over v1.0:
%                    - Auto-scaling eye vel y-axis per bout (not fixed ±500)
%                    - Shared FR y-axis across all bouts for comparison
%                    - FR clamped >= 0 (no negative smoothing artifacts)
%                    - xline for saccade markers (spans full panel height)
%                    - Spike raster ticks along top of each panel
%                    - min_bout_saccades default raised to 3
%                    - Better panel spacing with manual axes positioning
%                    - X tick labels only on bottom panel
%                    - Thicker eye velocity trace for visibility
% 2026-04-04  v1.0  Initial version.
