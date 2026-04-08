function [fig_handle, results, params] = plot_stationary_bout_snippets_ver_1_0(varargin)
% plot_stationary_bout_snippets_ver_1_0  Visualize stationary saccade bouts
% with eye velocity, firing rate, and saccade onset markers.
%
%   For a single neuron, finds the best platform-stationary periods
%   (ranked by total saccade count), and plots time-series snippets
%   showing eye velocity + smoothed firing rate + saccade onset markers.
%
%   Uses Chenhao's revised saccade onset times (saccade_ONSET_all_unit.mat)
%   which mark the true start of each saccade, not the velocity threshold
%   crossing used in Jeff's original saccade-edited.mat files.
%
%   INPUTS (name-value pairs via varargin):
%       neuron_name         - cell name (e.g., 'M402-2022-11-02-2-TT02_3')
%       neuron_index        - index into NPH_all_tfilelist (alternative)
%       tfilelist_path      - path to NPH all tfilelist.mat
%       chenhao_sacc_path   - path to saccade_ONSET_all_unit.mat
%       n_bouts             - how many best bouts to show (default: 3)
%       snippet_duration    - max seconds per panel (default: 15)
%       fr_bin_dt           - firing rate bin width in s (default: 0.01)
%       fr_smooth_win       - Gaussian smoothing window in bins (default: 15)
%       doSavePNG, doSaveFIG, doSaveEPS, doShow
%
%   OUTPUTS:
%       fig_handle  - figure handle ([] if insufficient data)
%       results     - struct with bout info, saccade counts, etc.
%       params      - struct of parameter values used
%
%   USAGE:
%       % By neuron name:
%       [fh, res, par] = plot_stationary_bout_snippets_ver_1_0(...
%           'neuron_name', 'M402-2022-11-02-2-TT02_3', 'doShow', 1);
%
%       % By index:
%       [fh, res, par] = plot_stationary_bout_snippets_ver_1_0(...
%           'neuron_index', 57, 'doShow', 1);
%
%       % Batch — top 5 cells from v1.3 best cells list:
%       top_cells = {'M402-2022-11-02-2-TT02_3', 'M284-2022-05-12-1-TT07_1', ...
%                    'M287-2022-08-08-TT08_1', 'M391-2022-12-07-1-TT01_1', ...
%                    'M247-2022-01-10-TT01_1'};
%       for ii = 1:numel(top_cells)
%           plot_stationary_bout_snippets_ver_1_0('neuron_name', top_cells{ii}, 'doShow', 1);
%       end
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
min_bout_saccades   = 2;        % minimum saccades in a bout to qualify

% Firing rate computation
fr_bin_dt           = 0.01;     % s
fr_smooth_win       = 15;       % bins for Gaussian smoothing

% Figure parameters
FontSize            = 20;
TickFontSize        = 16;
TitleFontSize       = 18;
LineWidth_EV        = 0.8;
LineWidth_FR        = 3;
LineWidth_sacc      = 2;
color_CCW_sacc      = [0.30 0.69 0.29];   % green — temporal / CCW
color_CW_sacc       = [0.90 0.45 0.15];   % orange — nasal / CW
color_FR            = [0.15 0.35 0.70];    % blue — firing rate
color_EV            = [0.40 0.40 0.40];    % grey — eye velocity
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
params.version              = '1.0';
params.date_run             = datestr(now, 'yyyy-mm-dd HH:MM:SS');
params.n_bouts              = n_bouts;
params.snippet_duration     = snippet_duration;
params.min_bout_saccades    = min_bout_saccades;
params.fr_bin_dt            = fr_bin_dt;
params.fr_smooth_win        = fr_smooth_win;
params.saccade_source       = 'Chenhao onset times (saccade_ONSET_all_unit.mat)';

%% ========================================================================
%  RESOLVE NEURON
%  ========================================================================
tfl = load(tfilelist_path, 'NPH_all_tfilelist');
NPH_all_tfilelist = tfl.NPH_all_tfilelist(:);

if ~isempty(neuron_index) && isempty(neuron_name)
    neuron_name_from_path = '';
    [~, neuron_name_from_path] = fileparts(NPH_all_tfilelist{neuron_index});
    neuron_name = neuron_name_from_path;
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
if verbose; fprintf('\n=== Stationary Bout Snippets ===\n'); end
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
            if verbose; fprintf('  Eye velocity: diffHdeg (deg/s)\n'); end
        elseif isfield(eye_fields, 'diffH')
            eye_vel_tsd = eye_fields.diffH;
            has_eye_data = true;
            if verbose; fprintf('  Eye velocity: diffH (pixels/s)\n'); end
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
ccw_per_period  = zeros(n_periods, 1);
cw_per_period   = zeros(n_periods, 1);

for iP = 1:n_periods
    nT_p = sum(temporal_stat > STtstart(iP) & temporal_stat < STtend(iP));
    nN_p = sum(nasal_stat > STtstart(iP) & nasal_stat < STtend(iP));
    sacc_per_period(iP) = nT_p + nN_p;
    ccw_per_period(iP) = nT_p;
    cw_per_period(iP)  = nN_p;
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
%  FIGURE
%  ========================================================================
fig_h = 280 * n_actual_bouts + 180;
if doShow
    fig_handle = figure('Position', [50 50 figWidth fig_h]);
else
    fig_handle = figure('Position', [50 50 figWidth fig_h], 'Visible', 'off');
end
fig_handle.UserData.created   = datestr(now, 'yyyy-mm-dd HH:MM:SS');
fig_handle.UserData.function  = 'plot_stationary_bout_snippets_ver_1_0';
fig_handle.UserData.neuron    = neuron_name;
fig_handle.UserData.saccade_source = 'Chenhao onset times';

% Ensure eye data vectors are columns
if has_eye_data
    ev_tvec = eye_vel_tsd.tvec(:);
    ev_data = eye_vel_tsd.data(:);
end

for iB = 1:n_actual_bouts
    pIdx = bout_indices(iB);
    t_win_start = STtstart(pIdx);
    t_win_end   = STtend(pIdx);
    period_dur  = t_win_end - t_win_start;

    % Get saccades in this period
    sacc_T_bout = temporal_stat(temporal_stat > t_win_start & temporal_stat < t_win_end);
    sacc_N_bout = nasal_stat(nasal_stat > t_win_start & nasal_stat < t_win_end);

    % If period longer than snippet_duration, center on densest saccade cluster
    if period_dur > snippet_duration
        all_sacc_in = sort([sacc_T_bout; sacc_N_bout]);
        if ~isempty(all_sacc_in)
            mid_sacc = median(all_sacc_in);
            t_win_start = max(STtstart(pIdx), mid_sacc - snippet_duration/2);
            t_win_end   = min(STtend(pIdx), t_win_start + snippet_duration);
        end
        % Re-filter saccades to trimmed window
        sacc_T_bout = sacc_T_bout(sacc_T_bout > t_win_start & sacc_T_bout < t_win_end);
        sacc_N_bout = sacc_N_bout(sacc_N_bout > t_win_start & sacc_N_bout < t_win_end);
    end

    ax = subplot(n_actual_bouts, 1, iB); hold(ax, 'on');

    % --- Eye velocity (left y-axis) ---
    yyaxis(ax, 'left');
    if has_eye_data
        ev_mask = ev_tvec >= t_win_start & ev_tvec <= t_win_end;
        if any(ev_mask)
            plot(ax, ev_tvec(ev_mask), ev_data(ev_mask), '-', ...
                'Color', color_EV, 'LineWidth', LineWidth_EV);
        end
    end
    ylabel(ax, 'Eye Vel (\circ/s)', 'FontSize', FontSize);
    ax.YAxis(1).Color = [0 0 0];
    ax.YAxis(1).FontSize = TickFontSize;

    % --- Saccade onset markers (vertical lines) ---
    % Draw these on left y-axis so they span the eye velocity range
    yl_left = ax.YLim;
    for iS = 1:length(sacc_T_bout)
        plot(ax, [sacc_T_bout(iS) sacc_T_bout(iS)], yl_left, '-', ...
            'Color', [color_CCW_sacc 0.65], 'LineWidth', LineWidth_sacc);
    end
    for iS = 1:length(sacc_N_bout)
        plot(ax, [sacc_N_bout(iS) sacc_N_bout(iS)], yl_left, '-', ...
            'Color', [color_CW_sacc 0.65], 'LineWidth', LineWidth_sacc);
    end

    % --- Firing rate (right y-axis) ---
    yyaxis(ax, 'right');
    fr_bin_edges = t_win_start:fr_bin_dt:t_win_end;
    if length(fr_bin_edges) > 2
        spk_in_win = myCell.t{1}(myCell.t{1} >= t_win_start & myCell.t{1} <= t_win_end);
        fr_bin_counts = histcounts(spk_in_win, fr_bin_edges);
        fr_inst_rate = fr_bin_counts / fr_bin_dt;
        fr_inst_smooth = smoothdata(fr_inst_rate, 'gaussian', fr_smooth_win);
        fr_bin_centers = fr_bin_edges(1:end-1) + fr_bin_dt/2;
        plot(ax, fr_bin_centers, fr_inst_smooth, '-', ...
            'Color', color_FR, 'LineWidth', LineWidth_FR);
    end
    ylabel(ax, 'FR (Hz)', 'FontSize', FontSize);
    ax.YAxis(2).Color = color_FR;
    ax.YAxis(2).FontSize = TickFontSize;

    % --- Formatting ---
    xlim(ax, [t_win_start t_win_end]);
    ax.FontSize = TickFontSize;
    box(ax, 'off');

    title_str = sprintf('Bout %d  —  %.1fs  |  %d CCW + %d CW saccades (Chenhao onsets)', ...
        iB, period_dur, length(sacc_T_bout), length(sacc_N_bout));
    title(ax, title_str, 'FontSize', TitleFontSize, 'FontWeight', 'bold');

    if iB == n_actual_bouts
        xlabel(ax, 'Time (s)', 'FontSize', FontSize);
    end
end

% --- Supertitle ---
sgtitle(sprintf('%s — Stationary Saccade Bouts', strrep(neuron_name, '_', '\_')), ...
    'FontSize', TitleFontSize + 4, 'FontWeight', 'bold');

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
if doSavePNG; print(fig_handle, [filepath_base '.png'], '-dpng', '-r600'); end
if doSaveFIG; savefig(fig_handle, [filepath_base '.fig']); end
if doSaveEPS; print(fig_handle, [filepath_base '.eps'], '-depsc', '-painters'); end

save(fullfile(out_folder, 'snippet_results.mat'), 'results', 'params', '-v7');

if verbose
    fprintf('\n  Saved to: %s\n', out_folder);
end

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
    % Fallback: match without session prefix
    for cc = 1:length(spike_cell_names)
        if contains(spike_cell_names{cc}, strrep(neuron_name, [SSN '-'], ''))
            cell_idx = cc; return;
        end
    end
end

%% Version History
% 2026-04-04  v1.0  Initial version. Standalone stationary bout snippet
%                    viewer. Uses Chenhao's saccade onset times from
%                    saccade_ONSET_all_unit.mat. Shows eye velocity +
%                    smoothed FR + CCW/CW saccade markers for the best
%                    stationary bouts per neuron.
