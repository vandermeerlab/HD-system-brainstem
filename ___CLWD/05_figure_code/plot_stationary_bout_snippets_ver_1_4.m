function [fig_handles, results, params] = plot_stationary_bout_snippets_ver_1_4(varargin)
% plot_stationary_bout_snippets_ver_1_4  Visualize stationary saccade bouts
%
%   For a single neuron, finds the best platform-stationary periods
%   (ranked by total saccade count) and plots ONE FIGURE PER BOUT with
%   4 vertically stacked sub-panels:
%     1. Eye position (deg)
%     2. Eye velocity (deg/s) with saccade onset markers
%     3. Smoothed firing rate (Hz) — computed via MakeQfromS pipeline
%     4. Spike train raster
%
%   Firing rate is computed on the FULL unrestricted spike train via
%   MakeQfromS with Gaussian smoothing, then restricted to each bout
%   window. This matches the canonical lab FR computation pattern.
%
%   For optotagged sessions, blue shaded bars mark laser stimulation.
%   Uses Chenhao's revised saccade onset times.
%
%   INPUTS (name-value pairs via varargin):
%       neuron_name         - cell name (e.g., 'M402-2022-11-02-2-TT02_3')
%       neuron_index        - index into NPH_all_tfilelist (alternative)
%       tfilelist_path      - path to NPH all tfilelist.mat
%       chenhao_sacc_path   - path to saccade_ONSET_all_unit.mat
%       n_bouts             - how many best bouts to show (default: 3)
%       snippet_duration    - max seconds per panel (default: 15)
%       min_bout_saccades   - minimum saccades per bout to qualify (default: 3)
%       fr_dt               - MakeQfromS bin width in s (default: 0.01)
%       fr_gausswin_sd      - Gaussian SD in s for MakeQfromS (default: 0.05)
%       ev_padding          - fraction of range to pad eye y-axes (default: 0.15)
%       doSavePNG, doSaveFIG, doSaveEPS, doShow
%
%   OUTPUTS:
%       fig_handles - cell array of figure handles (one per bout)
%       results     - struct with bout info, saccade counts, etc.
%       params      - struct of parameter values used
%
%   USAGE:
%       [fh, res, par] = plot_stationary_bout_snippets_ver_1_4(...
%           'neuron_name', 'M402-2022-11-02-2-TT02_3', 'doShow', 1);
%
%   DEPENDENCIES:
%       LoadSpikesQuick, SelectTS, MakeQfromS, restrict, tsd,
%       pushdir, popdir, HD_GetSSN, process_varargin,
%       SortBrainstemEventLabels3, load_chenhao_saccades,
%       save_fig_trio, find_cell_index, gausskernel, trim_histc

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

% Firing rate computation (MakeQfromS pipeline)
fr_dt               = 0.01;     % s — bin width
fr_gausswin_sd      = 0.05;     % s — Gaussian SD (50 ms, matches megaplot)

% Eye display
ev_padding          = 0.15;     % fraction of data range to pad y-axes

% Known optotagged sessions (for laser shading)
opto_sessions = {'M271-2021-08-30', 'M282-2022-02-10', 'M282-2022-02-15', ...
    'M391-2022-12-07-1', 'M391-2022-12-07-2', 'M402-2022-11-02-2', ...
    'M444-2023-08-18'};

% Figure parameters
FontSize            = 22;
TickFontSize        = 16;
TitleFontSize       = 20;
LabelFontSize       = 18;
LineWidth_EP        = 1.2;
LineWidth_EV        = 1.2;
LineWidth_FR        = 3.0;
LineWidth_sacc      = 2.0;
alpha_sacc          = 0.55;
color_CCW_sacc      = [0.30 0.69 0.29];   % green — temporal / CCW
color_CW_sacc       = [0.90 0.45 0.15];   % orange — nasal / CW
color_FR            = [0.85 0.15 0.15];    % red — firing rate
color_EP            = [0.10 0.10 0.10];    % near-black — eye position
color_EV            = [0.20 0.20 0.20];    % dark grey — eye velocity
color_spk           = [0.10 0.10 0.10];    % near-black — spike ticks
color_laser         = [0.30 0.50 0.90];    % blue — laser shading
alpha_laser         = 0.15;
figWidth            = 1600;
figHeight           = 750;

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
params.functionName         = 'plot_stationary_bout_snippets_ver_1_4';
params.dateRun              = datestr(now, 'yyyy-mm-dd HH:MM:SS');
params.saveRoot             = save_root;
params.nBouts               = n_bouts;
params.snippetDuration      = snippet_duration;
params.minBoutSaccades      = min_bout_saccades;
params.frDt                 = fr_dt;
params.frGausswinSd         = fr_gausswin_sd;
params.evPadding            = ev_padding;
params.saccadeSource        = 'Chenhao onset times (saccade_ONSET_all_unit.mat)';
params.frPipeline           = 'MakeQfromS (gauss smooth on full spike train, then restrict)';

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

params.neuronName = neuron_name;
if verbose; fprintf('\n=== Stationary Bout Snippets (v1.4) ===\n'); end
if verbose; fprintf('Neuron: %s\n', neuron_name); end

%% ========================================================================
%  LOAD SESSION DATA
%  ========================================================================
pushdir(sess_path);
SSN = HD_GetSSN;
params.SSN = SSN;
if verbose; fprintf('Session: %s\n', SSN); end

% --- Load Chenhao saccade onsets ---
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

% --- Compute FR via MakeQfromS on FULL spike train ---
cfg_Q = [];
cfg_Q.dt = fr_dt;
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = fr_gausswin_sd;
Q_full = MakeQfromS(cfg_Q, myCell);
FR_full = tsd(Q_full.tvec, Q_full.data ./ fr_dt);   % counts → Hz

if verbose; fprintf('  FR computed via MakeQfromS (dt=%.3f, gaussSD=%.3f)\n', fr_dt, fr_gausswin_sd); end

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
laser_on_times = []; laser_off_times = [];
if is_opto
    try
        [~, ~, laser_on_raw, laser_off_raw] = SortBrainstemEventLabels3;
        if ~(isscalar(laser_on_raw) && isnan(laser_on_raw))
            laser_on_times = laser_on_raw(:);
            if ~(isscalar(laser_off_raw) && isnan(laser_off_raw))
                laser_off_times = laser_off_raw(:);
            else
                laser_off_times = laser_on_times + 1.0;
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
    fig_handles = {}; results = struct(); return
end
n_actual_bouts = length(bout_indices);

%% ========================================================================
%  STORE RESULTS
%  ========================================================================
results = struct();
results.neuronName      = neuron_name;
results.SSN             = SSN;
results.nCCWStat        = length(temporal_stat);
results.nCWStat         = length(nasal_stat);
results.nStatPeriods    = n_periods;
results.statDurTotal    = sum(STtend - STtstart);
results.boutIndices     = bout_indices;
results.saccPerBout     = sacc_per_period(bout_indices);
results.saccadeSource   = 'Chenhao onset times';
results.isOpto          = is_opto;

%% ========================================================================
%  PRECOMPUTE DATA FOR EACH BOUT
%  ========================================================================
bout_data = cell(n_actual_bouts, 1);

if has_eye_vel
    ev_tvec = eye_vel_tsd.tvec(:); ev_data_all = eye_vel_tsd.data(:);
end
if has_eye_pos
    ep_tvec = eye_pos_tsd.tvec(:); ep_data_all = eye_pos_tsd.data(:);
end

for iB = 1:n_actual_bouts
    bd = struct();
    pIdx = bout_indices(iB);
    bd.tStartFull = STtstart(pIdx);
    bd.tEndFull   = STtend(pIdx);
    bd.periodDur  = bd.tEndFull - bd.tStartFull;

    % Get saccades in this period
    bd.saccT = temporal_stat(temporal_stat > bd.tStartFull & temporal_stat < bd.tEndFull);
    bd.saccN = nasal_stat(nasal_stat > bd.tStartFull & nasal_stat < bd.tEndFull);

    % Determine display window
    t_win_start = bd.tStartFull;
    t_win_end   = bd.tEndFull;
    if bd.periodDur > snippet_duration
        all_sacc_in = sort([bd.saccT; bd.saccN]);
        if ~isempty(all_sacc_in)
            mid_sacc = median(all_sacc_in);
            t_win_start = max(bd.tStartFull, mid_sacc - snippet_duration/2);
            t_win_end   = min(bd.tEndFull, t_win_start + snippet_duration);
        end
    end
    bd.tStart = t_win_start;
    bd.tEnd   = t_win_end;

    % Re-filter saccades to display window
    bd.saccTWin = bd.saccT(bd.saccT > t_win_start & bd.saccT < t_win_end);
    bd.saccNWin = bd.saccN(bd.saccN > t_win_start & bd.saccN < t_win_end);

    % Eye position in window
    bd.hasEP = false;
    if has_eye_pos
        ep_mask = ep_tvec >= t_win_start & ep_tvec <= t_win_end;
        if any(ep_mask)
            bd.epT = ep_tvec(ep_mask); bd.epD = ep_data_all(ep_mask);
            bd.hasEP = true;
        end
    end

    % Eye velocity in window
    bd.hasEV = false;
    if has_eye_vel
        ev_mask = ev_tvec >= t_win_start & ev_tvec <= t_win_end;
        if any(ev_mask)
            bd.evT = ev_tvec(ev_mask); bd.evD = ev_data_all(ev_mask);
            bd.hasEV = true;
        end
    end

    % Laser pulses in window
    bd.laserOnWin = []; bd.laserOffWin = [];
    if ~isempty(laser_on_times)
        laser_mask = laser_on_times >= t_win_start & laser_on_times <= t_win_end;
        if any(laser_mask)
            bd.laserOnWin  = laser_on_times(laser_mask);
            bd.laserOffWin = laser_off_times(laser_mask);
            bd.laserOffWin = min(bd.laserOffWin, t_win_end);
        end
    end

    % Firing rate — restrict precomputed MakeQfromS FR to window
    FR_bout = restrict(FR_full, t_win_start, t_win_end);
    bd.hasFR = false;
    if ~isempty(FR_bout.tvec)
        bd.frTvec   = FR_bout.tvec;
        bd.frData   = FR_bout.data;
        % Clamp negative smoothing artifacts
        bd.frData(bd.frData < 0) = 0;
        bd.hasFR    = true;
    end

    % Spike times in window
    spk_all = myCell.t{1};
    bd.spkTimes = spk_all(spk_all >= t_win_start & spk_all <= t_win_end);

    bout_data{iB} = bd;
end

% --- Shared y-axis limits across all bout figures ---
all_fr_max = 0; all_ep_min = Inf; all_ep_max = -Inf;
for iB = 1:n_actual_bouts
    bd = bout_data{iB};
    if bd.hasFR; all_fr_max = max(all_fr_max, max(bd.frData)); end
    if bd.hasEP
        all_ep_min = min(all_ep_min, min(bd.epD));
        all_ep_max = max(all_ep_max, max(bd.epD));
    end
end
fr_ylim_upper = ceil(all_fr_max / 10) * 10;
if fr_ylim_upper < 10; fr_ylim_upper = 10; end
ep_range = all_ep_max - all_ep_min;
if ep_range < 5; ep_range = 40; end
ep_ylim = [all_ep_min - ev_padding*ep_range, all_ep_max + ev_padding*ep_range];

%% ========================================================================
%  CREATE OUTPUT FOLDER
%  ========================================================================
date_str = datestr(now, 'yyyy-mm-dd');
run_num = 1;
out_folder = fullfile(save_root, sprintf('bout_snippets_%s_run%d', date_str, run_num));
while exist(out_folder, 'dir')
    run_num = run_num + 1;
    out_folder = fullfile(save_root, sprintf('bout_snippets_%s_run%d', date_str, run_num));
end
mkdir(out_folder);

%% ========================================================================
%  FIGURES — one per bout
%  ========================================================================
%  Layout (manual axes positioning):
%    Panel 1: Eye Position   — 22% of usable height
%    Panel 2: Eye Velocity   — 22% of usable height (saccade markers)
%    Panel 3: Firing Rate    — 36% of usable height
%    Panel 4: Spike Train    — 8%  of usable height
%    Inter-panel gaps        — 12% total

% Margins (normalized)
lm = 0.08;  rm = 0.03;  tm = 0.08;  bm = 0.08;
sub_gap = 0.015;

usable_h = 1 - tm - bm;
panel_w  = 1 - lm - rm;

% Panel height fractions (of usable height, after removing gaps)
frac_ep  = 0.22;
frac_ev  = 0.22;
frac_fr  = 0.36;
frac_spk = 0.08;
total_frac = frac_ep + frac_ev + frac_fr + frac_spk;
n_sub_gaps = 3;
avail_for_panels = usable_h - sub_gap * n_sub_gaps;

h_ep  = avail_for_panels * (frac_ep  / total_frac);
h_ev  = avail_for_panels * (frac_ev  / total_frac);
h_fr  = avail_for_panels * (frac_fr  / total_frac);
h_spk = avail_for_panels * (frac_spk / total_frac);

% Sub-panel y-positions (top to bottom: EP, EV, FR, SPK)
y_ep  = 1 - tm - h_ep;
y_ev  = y_ep  - sub_gap - h_ev;
y_fr  = y_ev  - sub_gap - h_fr;
y_spk = y_fr  - sub_gap - h_spk;

opto_str = '';
if is_opto; opto_str = '  [optotagged]'; end

fig_handles = cell(n_actual_bouts, 1);

for iB = 1:n_actual_bouts
    bd = bout_data{iB};
    x_lim = [bd.tStart bd.tEnd];

    % --- Create figure ---
    if doShow
        fh = figure('Position', [30 + 30*(iB-1) 50 figWidth figHeight]);
    else
        fh = figure('Position', [30 50 figWidth figHeight], 'Visible', 'off');
    end

    % Store metadata
    fh.UserData.created   = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fh.UserData.function  = params.functionName;
    fh.UserData.neuron    = neuron_name;
    fh.UserData.boutIndex = iB;
    fh.UserData.isOpto    = is_opto;
    fh.UserData.params    = params;

    all_axes = [];

    % ---- (1) Eye Position ----
    ax_ep = axes('Position', [lm y_ep panel_w h_ep]); hold(ax_ep, 'on'); %#ok<LAXES>
    if bd.hasEP
        plot(ax_ep, bd.epT, bd.epD, '-', 'Color', color_EP, 'LineWidth', LineWidth_EP);
        ylim(ax_ep, ep_ylim);
    end
    ylabel(ax_ep, 'Eye Pos (\circ)', 'FontSize', LabelFontSize);
    set(ax_ep, 'XTickLabel', [], 'FontSize', TickFontSize, 'Box', 'off');
    xlim(ax_ep, x_lim);
    title(ax_ep, sprintf('Bout %d  |  %.1fs  |  %d CCW + %d CW saccades', ...
        iB, bd.periodDur, length(bd.saccTWin), length(bd.saccNWin)), ...
        'FontSize', TitleFontSize, 'FontWeight', 'bold', 'Interpreter', 'none');
    all_axes(end+1) = ax_ep; %#ok<AGROW>

    % ---- (2) Eye Velocity + saccade markers ----
    ax_ev = axes('Position', [lm y_ev panel_w h_ev]); hold(ax_ev, 'on'); %#ok<LAXES>
    if bd.hasEV
        plot(ax_ev, bd.evT, bd.evD, '-', 'Color', color_EV, 'LineWidth', LineWidth_EV);
        ev_range_bout = max(bd.evD) - min(bd.evD);
        if ev_range_bout < 10; ev_range_bout = 100; end
        ev_lo = min(bd.evD) - ev_padding * ev_range_bout;
        ev_hi = max(bd.evD) + ev_padding * ev_range_bout;
        ylim(ax_ev, [ev_lo ev_hi]);
    end
    % Saccade onset markers (CCW = green, CW = orange)
    for iS = 1:length(bd.saccTWin)
        xline(ax_ev, bd.saccTWin(iS), '-', 'Color', [color_CCW_sacc alpha_sacc], ...
            'LineWidth', LineWidth_sacc);
    end
    for iS = 1:length(bd.saccNWin)
        xline(ax_ev, bd.saccNWin(iS), '-', 'Color', [color_CW_sacc alpha_sacc], ...
            'LineWidth', LineWidth_sacc);
    end
    ylabel(ax_ev, 'Eye Vel (\circ/s)', 'FontSize', LabelFontSize);
    set(ax_ev, 'XTickLabel', [], 'FontSize', TickFontSize, 'Box', 'off');
    xlim(ax_ev, x_lim);
    all_axes(end+1) = ax_ev; %#ok<AGROW>

    % ---- (3) Firing Rate ----
    ax_fr = axes('Position', [lm y_fr panel_w h_fr]); hold(ax_fr, 'on'); %#ok<LAXES>
    if bd.hasFR
        plot(ax_fr, bd.frTvec, bd.frData, '-', 'Color', color_FR, 'LineWidth', LineWidth_FR);
        ylim(ax_fr, [0 fr_ylim_upper]);
    end
    % Saccade onset markers on FR panel too (subtle, thinner)
    for iS = 1:length(bd.saccTWin)
        xline(ax_fr, bd.saccTWin(iS), '-', 'Color', [color_CCW_sacc 0.25], ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    for iS = 1:length(bd.saccNWin)
        xline(ax_fr, bd.saccNWin(iS), '-', 'Color', [color_CW_sacc 0.25], ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    ylabel(ax_fr, 'FR (Hz)', 'FontSize', LabelFontSize);
    set(ax_fr, 'XTickLabel', [], 'FontSize', TickFontSize, 'Box', 'off');
    xlim(ax_fr, x_lim);
    all_axes(end+1) = ax_fr; %#ok<AGROW>

    % ---- (4) Spike Train ----
    ax_spk = axes('Position', [lm y_spk panel_w h_spk]); hold(ax_spk, 'on'); %#ok<LAXES>
    if ~isempty(bd.spkTimes)
        spk_x = repmat(bd.spkTimes(:)', 2, 1);
        spk_y = repmat([0; 1], 1, length(bd.spkTimes));
        plot(ax_spk, spk_x, spk_y, '-', 'Color', color_spk, 'LineWidth', 0.5);
    end
    ylim(ax_spk, [0 1]);
    set(ax_spk, 'YTick', [], 'FontSize', TickFontSize, 'Box', 'off');
    xlim(ax_spk, x_lim);
    xlabel(ax_spk, 'Time (s)', 'FontSize', FontSize);
    ylabel(ax_spk, 'Spikes', 'FontSize', LabelFontSize - 4);
    all_axes(end+1) = ax_spk; %#ok<AGROW>

    % ---- Laser shading (all sub-panels) ----
    if ~isempty(bd.laserOnWin)
        for iL = 1:length(bd.laserOnWin)
            l_on  = bd.laserOnWin(iL);
            l_off = bd.laserOffWin(iL);
            for iAx = 1:length(all_axes)
                yl_ax = ylim(all_axes(iAx));
                nan_check = [l_on l_off yl_ax(1) yl_ax(2)];
                if any(isnan(nan_check)); continue; end
                fill(all_axes(iAx), [l_on l_off l_off l_on], ...
                    [yl_ax(1) yl_ax(1) yl_ax(2) yl_ax(2)], ...
                    color_laser, 'EdgeColor', 'none', 'FaceAlpha', alpha_laser, ...
                    'HandleVisibility', 'off');
            end
        end
    end

    % --- Supertitle ---
    sgtitle(sprintf('%s  |  Bout %d of %d%s', ...
        neuron_name, iB, n_actual_bouts, opto_str), ...
        'FontSize', FontSize, 'FontWeight', 'bold', 'Interpreter', 'none');

    % --- Save this bout figure ---
    safe_name = strrep(neuron_name, '-', '_');
    filepath_base = fullfile(out_folder, sprintf('snippet_%s_bout%d', safe_name, iB));
    save_fig_trio(fh, filepath_base, ...
        'doSavePNG', doSavePNG, 'doSaveFIG', doSaveFIG, 'doSaveEPS', doSaveEPS, ...
        'doShow', doShow);

    if doShow
        fig_handles{iB} = fh;
    else
        fig_handles{iB} = [];
    end

    if verbose; fprintf('  Saved bout %d/%d\n', iB, n_actual_bouts); end
end

%% ========================================================================
%  SAVE RESULTS
%  ========================================================================
save(fullfile(out_folder, 'snippet_results.mat'), 'results', 'params', '-v7');

if verbose; fprintf('\n  Output folder: %s\n', out_folder); end
if verbose; fprintf('  Done.\n'); end

end

%% Version History
% 2026-04-06  v1.4  Major rewrite:
%                    - ONE FIGURE PER BOUT (was all bouts stacked)
%                    - FR via MakeQfromS pipeline (was histcounts +
%                      smoothdata). Computes on full spike train, then
%                      restricts to each bout window. Matches canonical
%                      lab FR computation.
%                    - Subtle saccade markers on FR panel (thin, low alpha)
%                    - Vectorized spike raster drawing (was per-spike loop)
%                    - Bigger default fonts (22/20/18/16)
%                    - Wider figure (1600px) with better vertical spacing
%                    - params camelCase throughout
%                    - Each bout saved as separate file:
%                      snippet_{neuron}_bout{N}.png/fig/eps
%                    - NaN guard on laser fill()
% 2026-04-05  v1.3  4 stacked sub-panels per bout, blue laser shading,
%                    red FR trace, HandleVisibility off, Interpreter none
% 2026-04-05  v1.2  Shared helpers refactor
% 2026-04-04  v1.1  Auto-scaling, shared FR axis, xline markers
% 2026-04-04  v1.0  Initial version
