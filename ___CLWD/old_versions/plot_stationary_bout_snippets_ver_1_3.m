function [fig_handle, results, params] = plot_stationary_bout_snippets_ver_1_3(varargin)
% plot_stationary_bout_snippets_ver_1_3  Visualize stationary saccade bouts
% with stacked time-aligned panels: eye position, eye velocity, firing
% rate, and spike train.
%
%   For a single neuron, finds the best platform-stationary periods
%   (ranked by total saccade count), and plots multi-panel time-series
%   snippets. Each bout gets 4 vertically stacked sub-panels:
%     1. Eye position (deg)
%     2. Eye velocity (deg/s) with saccade onset markers
%     3. Smoothed firing rate (Hz)
%     4. Spike train raster
%   For optotagged sessions, blue shaded bars mark laser stimulation periods.
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
%       ev_padding          - fraction of range to pad eye y-axes (default: 0.15)
%       doSavePNG, doSaveFIG, doSaveEPS, doShow
%
%   OUTPUTS:
%       fig_handle  - figure handle ([] if insufficient data)
%       results     - struct with bout info, saccade counts, etc.
%       params      - struct of parameter values used
%
%   USAGE:
%       [fh, res, par] = plot_stationary_bout_snippets_ver_1_3(...
%           'neuron_name', 'M402-2022-11-02-2-TT02_3', 'doShow', 1);
%
%   DEPENDENCIES:
%       LoadSpikesQuick, SelectTS, restrict, pushdir, popdir,
%       HD_GetSSN, process_varargin, SortBrainstemEventLabels3,
%       load_chenhao_saccades, save_fig_trio, find_cell_index

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

% Eye display
ev_padding          = 0.15;     % fraction of data range to pad y-axes

% Known optotagged sessions (for laser shading)
opto_sessions = {'M271-2021-08-30', 'M282-2022-02-10', 'M282-2022-02-15', ...
    'M391-2022-12-07-1', 'M391-2022-12-07-2', 'M402-2022-11-02-2', ...
    'M444-2023-08-18'};

% Figure parameters
FontSize            = 18;
TickFontSize        = 14;
TitleFontSize       = 16;
LabelFontSize       = 16;
LineWidth_EP        = 1.0;
LineWidth_EV        = 1.0;
LineWidth_FR        = 2.5;
LineWidth_sacc      = 1.8;
alpha_sacc          = 0.50;
color_CCW_sacc      = [0.30 0.69 0.29];   % green — temporal / CCW
color_CW_sacc       = [0.90 0.45 0.15];   % orange — nasal / CW
color_FR            = [0.85 0.15 0.15];    % red — firing rate (matches reference)
color_EP            = [0.10 0.10 0.10];    % black — eye position
color_EV            = [0.20 0.20 0.20];    % dark grey — eye velocity
color_spk           = [0.10 0.10 0.10];    % black — spike ticks
color_laser         = [0.30 0.50 0.90];    % blue — laser shading
alpha_laser         = 0.15;
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
params.functionName        = 'plot_stationary_bout_snippets_ver_1_3';
params.version              = '1.3';
params.dateRun             = datestr(now, 'yyyy-mm-dd HH:MM:SS');
params.saveRoot             = save_root;
params.n_bouts              = n_bouts;
params.snippet_duration     = snippet_duration;
params.min_bout_saccades    = min_bout_saccades;
params.fr_bin_dt            = fr_bin_dt;
params.fr_smooth_win        = fr_smooth_win;
params.ev_padding           = ev_padding;
params.saccadeSource       = 'Chenhao onset times (saccade_ONSET_all_unit.mat)';

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
if verbose; fprintf('\n=== Stationary Bout Snippets (v1.3) ===\n'); end
if verbose; fprintf('Neuron: %s\n', neuron_name); end

%% ========================================================================
%  LOAD SESSION DATA
%  ========================================================================
pushdir(sess_path);
SSN = HD_GetSSN;
params.SSN = SSN;
if verbose; fprintf('Session: %s\n', SSN); end

% --- Load Chenhao saccade onsets for this session ---
[temporal_all, nasal_all] = load_chenhao_saccades(SSN, ...
    'chenhao_sacc_path', chenhao_sacc_path);

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

% --- Load eye data (velocity + position) ---
sacc_file = strcat(SSN, '-saccades-edited.mat');
has_eye_vel = false; has_eye_pos = false;
eye_vel_tsd = []; eye_pos_tsd = [];
if exist(sacc_file, 'file')
    try
        eye_fields = load(sacc_file, 'diffHdeg', 'diffH', 'tsdH');
        if isfield(eye_fields, 'diffHdeg')
            eye_vel_tsd = eye_fields.diffHdeg; has_eye_vel = true;
        elseif isfield(eye_fields, 'diffH')
            eye_vel_tsd = eye_fields.diffH; has_eye_vel = true;
        end
        if isfield(eye_fields, 'tsdH')
            eye_pos_tsd = eye_fields.tsdH; has_eye_pos = true;
        end
    catch ME_eye
        if verbose; fprintf('  Could not load eye data: %s\n', ME_eye.message); end
    end
end
if verbose
    fprintf('  Eye data: velocity=%d, position=%d\n', has_eye_vel, has_eye_pos);
end

% --- Load laser times (optotagged sessions) ---
is_opto = any(strcmp(SSN, opto_sessions));
laser_on_times = [];
if is_opto
    try
        [~, ~, laser_on_raw, laser_off_raw] = SortBrainstemEventLabels3;
        if ~(isscalar(laser_on_raw) && isnan(laser_on_raw))
            laser_on_times = laser_on_raw(:);
            if ~(isscalar(laser_off_raw) && isnan(laser_off_raw))
                laser_off_times = laser_off_raw(:);
            else
                laser_off_times = laser_on_times + 1.0;  % assume 1s pulse
            end
        end
    catch ME_laser
        if verbose; fprintf('  Laser detection failed: %s\n', ME_laser.message); end
    end
end
if verbose && is_opto
    if ~isempty(laser_on_times)
        fprintf('  Laser: %d pulses detected\n', length(laser_on_times));
    else
        fprintf('  Laser: opto session but no laser times found\n');
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
    fig_handle = []; results = struct(); popdir; return
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
results.is_opto         = is_opto;

%% ========================================================================
%  PRECOMPUTE DATA FOR EACH BOUT
%  ========================================================================
bout_data = cell(n_actual_bouts, 1);

if has_eye_vel
    ev_tvec = eye_vel_tsd.tvec(:); ev_data = eye_vel_tsd.data(:);
end
if has_eye_pos
    ep_tvec = eye_pos_tsd.tvec(:); ep_data = eye_pos_tsd.data(:);
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

    % Eye position in window
    bd.has_ep = false;
    if has_eye_pos
        ep_mask = ep_tvec >= t_win_start & ep_tvec <= t_win_end;
        if any(ep_mask)
            bd.ep_t = ep_tvec(ep_mask); bd.ep_d = ep_data(ep_mask);
            bd.has_ep = true;
        end
    end

    % Eye velocity in window
    bd.has_ev = false;
    if has_eye_vel
        ev_mask = ev_tvec >= t_win_start & ev_tvec <= t_win_end;
        if any(ev_mask)
            bd.ev_t = ev_tvec(ev_mask); bd.ev_d = ev_data(ev_mask);
            bd.has_ev = true;
        end
    end

    % Laser pulses in window
    bd.laser_on_win = []; bd.laser_off_win = [];
    if ~isempty(laser_on_times)
        laser_mask = laser_on_times >= t_win_start & laser_on_times <= t_win_end;
        if any(laser_mask)
            bd.laser_on_win  = laser_on_times(laser_mask);
            bd.laser_off_win = laser_off_times(laser_mask);
            % Clip to window boundaries
            bd.laser_off_win = min(bd.laser_off_win, t_win_end);
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
        bd.fr_smooth(bd.fr_smooth < 0) = 0;
        bd.fr_centers = fr_bin_edges(1:end-1) + fr_bin_dt/2;
        bd.spk_times = spk_in_win;
        bd.has_fr = true;
    end

    bout_data{iB} = bd;
end

% --- Shared y-axis limits across bouts ---
all_fr_max = 0; all_ep_min = Inf; all_ep_max = -Inf;
for iB = 1:n_actual_bouts
    bd = bout_data{iB};
    if bd.has_fr; all_fr_max = max(all_fr_max, max(bd.fr_smooth)); end
    if bd.has_ep
        all_ep_min = min(all_ep_min, min(bd.ep_d));
        all_ep_max = max(all_ep_max, max(bd.ep_d));
    end
end
fr_ylim_upper = ceil(all_fr_max / 10) * 10;
if fr_ylim_upper < 10; fr_ylim_upper = 10; end
ep_range = all_ep_max - all_ep_min;
if ep_range < 5; ep_range = 40; end
ep_ylim = [all_ep_min - ev_padding*ep_range, all_ep_max + ev_padding*ep_range];

%% ========================================================================
%  FIGURE — stacked multi-panel layout
%  ========================================================================
%  Each bout has 4 sub-panels:
%    [1] Eye Position   — 22% of bout height
%    [2] Eye Velocity   — 22% of bout height (with saccade markers)
%    [3] Firing Rate    — 28% of bout height
%    [4] Spike Train    —  8% of bout height
%  Internal gap between sub-panels: 2% of bout height
%  Gap between bouts: 4% of total figure height

fig_h = 550 * n_actual_bouts + 120;
if doShow
    fig_handle = figure('Position', [30 30 figWidth fig_h]);
else
    fig_handle = figure('Position', [30 30 figWidth fig_h], 'Visible', 'off');
end
fig_handle.UserData.created        = datestr(now, 'yyyy-mm-dd HH:MM:SS');
fig_handle.UserData.function       = 'plot_stationary_bout_snippets_ver_1_3';
fig_handle.UserData.neuron         = neuron_name;
fig_handle.UserData.saccade_source = 'Chenhao onset times';
fig_handle.UserData.is_opto        = is_opto;

% Layout constants (normalized figure coordinates)
lm = 0.08; rm = 0.04; tm = 0.05; bm = 0.05;
bout_gap = 0.03;   % gap between bouts
sub_gap  = 0.008;  % gap between sub-panels within a bout

usable_h = 1 - tm - bm;
total_bout_gaps = bout_gap * (n_actual_bouts - 1);
bout_block_h = (usable_h - total_bout_gaps) / n_actual_bouts;
panel_w = 1 - lm - rm;

% Sub-panel proportions within each bout (must sum to ~1 after gaps)
frac_ep  = 0.24;   % eye position
frac_ev  = 0.24;   % eye velocity
frac_fr  = 0.32;   % firing rate
frac_spk = 0.08;   % spike train
total_frac = frac_ep + frac_ev + frac_fr + frac_spk;
n_sub_gaps = 3;     % gaps between 4 panels
avail_for_panels = bout_block_h - sub_gap * n_sub_gaps;

h_ep  = avail_for_panels * (frac_ep / total_frac);
h_ev  = avail_for_panels * (frac_ev / total_frac);
h_fr  = avail_for_panels * (frac_fr / total_frac);
h_spk = avail_for_panels * (frac_spk / total_frac);

for iB = 1:n_actual_bouts
    bd = bout_data{iB};
    x_lim = [bd.t_start bd.t_end];

    % Top of this bout block
    bout_top = 1 - tm - (iB - 1) * (bout_block_h + bout_gap);

    % Sub-panel y-positions (top to bottom: EP, EV, FR, SPK)
    y_ep  = bout_top - h_ep;
    y_ev  = y_ep - sub_gap - h_ev;
    y_fr  = y_ev - sub_gap - h_fr;
    y_spk = y_fr - sub_gap - h_spk;

    % Axes handles for laser shading
    all_axes = [];

    % ---- (1) Eye Position ----
    ax_ep = axes('Position', [lm y_ep panel_w h_ep]); hold(ax_ep, 'on'); %#ok<LAXES>
    if bd.has_ep
        plot(ax_ep, bd.ep_t, bd.ep_d, '-', 'Color', color_EP, 'LineWidth', LineWidth_EP);
        ylim(ax_ep, ep_ylim);
    end
    ylabel(ax_ep, 'Eye Pos (\circ)', 'FontSize', LabelFontSize);
    set(ax_ep, 'XTickLabel', [], 'FontSize', TickFontSize, 'Box', 'off');
    xlim(ax_ep, x_lim);
    title(ax_ep, sprintf('Bout %d  —  %.1fs  |  %d CCW + %d CW saccades', ...
        iB, bd.period_dur, length(bd.sacc_T_win), length(bd.sacc_N_win)), ...
        'FontSize', TitleFontSize, 'FontWeight', 'bold', 'Interpreter', 'none');
    all_axes(end+1) = ax_ep; %#ok<AGROW>

    % ---- (2) Eye Velocity + saccade markers ----
    ax_ev = axes('Position', [lm y_ev panel_w h_ev]); hold(ax_ev, 'on'); %#ok<LAXES>
    if bd.has_ev
        plot(ax_ev, bd.ev_t, bd.ev_d, '-', 'Color', color_EV, 'LineWidth', LineWidth_EV);
        ev_range_bout = max(bd.ev_d) - min(bd.ev_d);
        if ev_range_bout < 10; ev_range_bout = 100; end
        ev_lo = min(bd.ev_d) - ev_padding * ev_range_bout;
        ev_hi = max(bd.ev_d) + ev_padding * ev_range_bout;
        ylim(ax_ev, [ev_lo ev_hi]);
    end
    % Saccade onset markers
    for iS = 1:length(bd.sacc_T_win)
        xline(ax_ev, bd.sacc_T_win(iS), '-', 'Color', [color_CCW_sacc alpha_sacc], ...
            'LineWidth', LineWidth_sacc);
    end
    for iS = 1:length(bd.sacc_N_win)
        xline(ax_ev, bd.sacc_N_win(iS), '-', 'Color', [color_CW_sacc alpha_sacc], ...
            'LineWidth', LineWidth_sacc);
    end
    ylabel(ax_ev, 'Eye Vel (\circ/s)', 'FontSize', LabelFontSize);
    set(ax_ev, 'XTickLabel', [], 'FontSize', TickFontSize, 'Box', 'off');
    xlim(ax_ev, x_lim);
    all_axes(end+1) = ax_ev; %#ok<AGROW>

    % ---- (3) Firing Rate ----
    ax_fr = axes('Position', [lm y_fr panel_w h_fr]); hold(ax_fr, 'on'); %#ok<LAXES>
    if bd.has_fr
        plot(ax_fr, bd.fr_centers, bd.fr_smooth, '-', 'Color', color_FR, 'LineWidth', LineWidth_FR);
        ylim(ax_fr, [0 fr_ylim_upper]);
    end
    ylabel(ax_fr, 'FR (Hz)', 'FontSize', LabelFontSize);
    set(ax_fr, 'XTickLabel', [], 'FontSize', TickFontSize, 'Box', 'off');
    xlim(ax_fr, x_lim);
    all_axes(end+1) = ax_fr; %#ok<AGROW>

    % ---- (4) Spike Train ----
    ax_spk = axes('Position', [lm y_spk panel_w h_spk]); hold(ax_spk, 'on'); %#ok<LAXES>
    if bd.has_fr && ~isempty(bd.spk_times)
        for iSpk = 1:length(bd.spk_times)
            plot(ax_spk, [bd.spk_times(iSpk) bd.spk_times(iSpk)], [0 1], '-', ...
                'Color', color_spk, 'LineWidth', 0.5);
        end
    end
    ylim(ax_spk, [0 1]);
    set(ax_spk, 'YTick', [], 'FontSize', TickFontSize, 'Box', 'off');
    xlim(ax_spk, x_lim);
    if iB == n_actual_bouts
        xlabel(ax_spk, 'Time (s)', 'FontSize', FontSize);
    else
        set(ax_spk, 'XTickLabel', []);
    end
    ylabel(ax_spk, 'Spikes', 'FontSize', LabelFontSize - 2);
    all_axes(end+1) = ax_spk; %#ok<AGROW>

    % ---- Laser shading (all sub-panels) ----
    if ~isempty(bd.laser_on_win)
        for iL = 1:length(bd.laser_on_win)
            l_on  = bd.laser_on_win(iL);
            l_off = bd.laser_off_win(iL);
            for iAx = 1:length(all_axes)
                yl_ax = ylim(all_axes(iAx));
                fill(all_axes(iAx), [l_on l_off l_off l_on], ...
                    [yl_ax(1) yl_ax(1) yl_ax(2) yl_ax(2)], ...
                    color_laser, 'EdgeColor', 'none', 'FaceAlpha', alpha_laser, ...
                    'HandleVisibility', 'off');
            end
        end
    end
end

% --- Supertitle ---
opto_str = '';
if is_opto; opto_str = '  [optotagged]'; end
sgtitle(sprintf('%s — Stationary Saccade Bouts (Chenhao onsets)%s', ...
    neuron_name, opto_str), 'FontSize', TitleFontSize + 4, 'FontWeight', 'bold', ...
    'Interpreter', 'none');

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
save_fig_trio(fig_handle, filepath_base, ...
    'doSavePNG', doSavePNG, 'doSaveFIG', doSaveFIG, 'doSaveEPS', doSaveEPS, 'doShow', doShow);

save(fullfile(out_folder, 'snippet_results.mat'), 'results', 'params', '-v7');

if verbose; fprintf('\n  Saved to: %s\n', out_folder); end
if verbose; fprintf('  Done.\n'); end

end

%% Version History
% 2026-04-05  v1.3  Major layout overhaul:
%                    - 4 stacked sub-panels per bout (eye position, eye
%                      velocity, firing rate, spike train) replacing the
%                      single yyaxis panel from v1.2
%                    - Eye position panel (tsdH) — new, shared y-axis
%                      across bouts
%                    - Saccade onset markers on eye velocity panel only
%                    - Dedicated spike train panel (black ticks, no y-axis)
%                    - Blue shaded bars for laser stimulation periods on
%                      optotagged sessions (all 4 sub-panels)
%                    - FR trace color changed to red (matches reference fig)
%                    - Added LabelFontSize parameter for y-axis labels
%                    - HandleVisibility off on laser shading fill()
%                    - Interpreter none on all title/sgtitle calls
% 2026-04-05  v1.2  Refactor for shared utilities:
%                    - Chenhao saccade loading via load_chenhao_saccades.m
%                    - Figure saving via save_fig_trio.m (shared utility)
%                    - Cell matching via find_cell_index.m (shared utility)
%                    - Removed duplicated local helper functions
% 2026-04-04  v1.1  Improvements over v1.0:
%                    - Auto-scaling eye vel y-axis per bout (not fixed +/-500)
%                    - Shared FR y-axis across all bouts for comparison
%                    - FR clamped >= 0 (no negative smoothing artifacts)
%                    - xline for saccade markers (spans full panel height)
%                    - Spike raster ticks along top of each panel
%                    - min_bout_saccades default raised to 3
%                    - Better panel spacing with manual axes positioning
%                    - X tick labels only on bottom panel
%                    - Thicker eye velocity trace for visibility
% 2026-04-04  v1.0  Initial version.
