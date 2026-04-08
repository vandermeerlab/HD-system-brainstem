function [shuffle_out, params] = compute_optotag_eye_shuffle(ep_tvec, ep_data_deg, laser_on, session_duration, cfg_in, varargin)
% compute_optotag_eye_shuffle  Shuffle test for laser-evoked eye displacement
%
%   Tests whether laser stimulation of optotagged NPH->DTN neurons evokes
%   significant eye movements. Generates a null distribution of |eye
%   displacement| by placing laser onsets at random times, then compares
%   the real laser-aligned |displacement| against the null.
%
%   THREE SHUFFLE METHODS (cfg.shuffle_method):
%     'unrestricted'        — circular shift across entire session (v1.0)
%     'random_stationary'   — randomly place onsets within stationary
%                             periods (DEFAULT, recommended)
%     'circular_stationary' — circular shift within concatenated stationary
%                             epoch timeline, then map back to real time
%
%   For 'random_stationary' and 'circular_stationary', the stationary
%   period intervals must be provided via cfg.stat_intervals = [nInt x 2]
%   matrix of [start, end] times.
%
%   INPUTS:
%       ep_tvec          - [nSamples x 1] eye position time vector (s,
%                          session-relative, from tsdH.tvec)
%       ep_data_deg      - [nSamples x 1] eye position in degrees
%                          (already multiplied by calibration factor k)
%       laser_on         - [nTrials x 1] laser onset times (s, session-relative)
%       session_duration - scalar, total session duration in seconds
%       cfg_in           - config structure (all fields optional)
%       varargin         - name-value pairs overriding any cfg field
%
%   CFG FIELDS:
%       .nShuffles       - number of shuffles (default: 500)
%       .window          - alignment window in s (default: [-1 2])
%       .dt              - bin size for interpolated eye grid (default: 0.005)
%       .laser_dur       - laser pulse duration in s (default: 1.0)
%       .baseline        - [start end] for baseline subtraction (default: [-1 -0.2])
%       .smoothEye       - smooth mean displacement traces (default: 1)
%       .smooth_win      - Gaussian smoothing window in bins (default: 5)
%       .alpha_001       - threshold for p < .001 (default: 0.001)
%       .alpha_05        - threshold for p < .05 (default: 0.05)
%       .shuffle_method  - 'unrestricted', 'random_stationary', or
%                          'circular_stationary' (default: 'random_stationary')
%       .stat_intervals  - [nIntervals x 2] matrix of [start, end] times
%                          for stationary periods (required for stationary methods)
%       .verbose         - print progress (default: 1)
%
%   OUTPUTS:
%       shuffle_out     - struct with fields:
%           .eye_grid, .disp_real, .disp_real_smooth
%           .shuff_matrix, .shuff_mean, .shuff_sem
%           .pval_per_bin, .sig_001, .sig_05, .laser_bins
%           .peak_disp_real, .peak_disp_shuff, .percentile_overall
%           .n_valid_trials, .shuffle_method, .stat_duration_total
%       params          - parameter struct
%
%   METHODOLOGICAL NOTE (Claude-written, flagged to Jeff):
%       This function was written by Claude. The stationary-period
%       restriction ensures the null distribution reflects chance eye
%       displacement during platform-stationary epochs (when laser was
%       actually delivered), not VOR-driven eye movements.
%
%   DEPENDENCIES:
%       ProcessConfig2, process_varargin

%% Config defaults
cfg_def = [];
cfg_def.nShuffles       = 500;
cfg_def.window          = [-1 2];
cfg_def.dt              = 0.005;
cfg_def.laser_dur       = 1.0;
cfg_def.baseline        = [-1 -0.2];
cfg_def.smoothEye       = 1;
cfg_def.smooth_win      = 5;
cfg_def.alpha_001       = 0.001;
cfg_def.alpha_05        = 0.05;
cfg_def.test_window     = [0 0.2];   % window for significance testing (s from laser onset)
cfg_def.shuffle_method  = 'random_stationary';
cfg_def.stat_intervals  = [];   % [nInt x 2] matrix of [start, end]
cfg_def.verbose         = 1;

if nargin < 5 || isempty(cfg_in); cfg_in = []; end
cfg = ProcessConfig2(cfg_def, cfg_in);

% Manual extraction loop for struct/cell params in varargin
for iV = 1:2:length(varargin)
    if ischar(varargin{iV})
        switch varargin{iV}
            case 'stat_intervals'; cfg.stat_intervals = varargin{iV+1};
        end
    end
end
process_varargin(varargin);

%% Validate
ep_tvec     = ep_tvec(:);
ep_data_deg = ep_data_deg(:);
laser_on    = laser_on(:);

use_stationary = ismember(cfg.shuffle_method, {'random_stationary', 'circular_stationary'});
if use_stationary && isempty(cfg.stat_intervals)
    error('compute_optotag_eye_shuffle:noStatIntervals', ...
        'shuffle_method=''%s'' requires cfg.stat_intervals [nInt x 2].', cfg.shuffle_method);
end

%% Precompute stationary period info
stat_duration_total = 0;
if use_stationary
    stat_iv = cfg.stat_intervals;
    stat_durations = stat_iv(:,2) - stat_iv(:,1);
    stat_duration_total = sum(stat_durations);
    stat_cumsum = cumsum(stat_durations);  % for circular method
end

%% Store params
params                  = [];
params.function_name    = 'compute_optotag_eye_shuffle';
params.version          = '1.2';
params.date_run         = datestr(now, 'yyyy-mm-dd HH:MM:SS');
params.nShuffles        = cfg.nShuffles;
params.window           = cfg.window;
params.dt               = cfg.dt;
params.laser_dur        = cfg.laser_dur;
params.baseline         = cfg.baseline;
params.smoothEye        = cfg.smoothEye;
params.smooth_win       = cfg.smooth_win;
params.alpha_001        = cfg.alpha_001;
params.alpha_05         = cfg.alpha_05;
params.test_window      = cfg.test_window;
params.shuffle_method   = cfg.shuffle_method;
params.nTrials          = length(laser_on);
params.session_duration = session_duration;
params.nSamples_eye     = length(ep_tvec);
params.stat_duration_total = stat_duration_total;
if use_stationary
    params.n_stat_intervals = size(stat_iv, 1);
end

%% Set up interpolation grid
eye_grid = cfg.window(1) : cfg.dt : cfg.window(2);
eye_grid = eye_grid(1:end-1) + cfg.dt / 2;
nBins    = length(eye_grid);

%% Masks
% test_bins: used for significance testing (default [0, 0.2] — early transient)
% laser_bins kept for backward compat in output struct (full laser window)
test_bins  = (eye_grid >= cfg.test_window(1)) & (eye_grid < cfg.test_window(2));
laser_bins = (eye_grid >= 0) & (eye_grid < cfg.laser_dur);
bl_mask    = (eye_grid >= cfg.baseline(1)) & (eye_grid <= cfg.baseline(2));

nTrials = length(laser_on);

%% ================================================================
%  Compute REAL mean |displacement|
%  ================================================================
[disp_real, n_valid_real] = compute_mean_abs_displacement( ...
    ep_tvec, ep_data_deg, laser_on, eye_grid, bl_mask, cfg.window);

if cfg.smoothEye && cfg.smooth_win > 1
    disp_real_smooth = smoothdata(disp_real, 'gaussian', cfg.smooth_win);
else
    disp_real_smooth = disp_real;
end

if cfg.verbose
    fprintf('compute_optotag_eye_shuffle [%s]: %d valid trials (of %d)\n', ...
        cfg.shuffle_method, n_valid_real, nTrials);
    if use_stationary
        fprintf('  Stationary time: %.1f s across %d intervals\n', ...
            stat_duration_total, size(stat_iv, 1));
    end
end

%% ================================================================
%  Shuffle loop
%  ================================================================
shuff_matrix = NaN(cfg.nShuffles, nBins);

if cfg.verbose
    fprintf('  Running %d shuffles (%s)...\n', cfg.nShuffles, cfg.shuffle_method);
end

for iShuff = 1:cfg.nShuffles

    switch cfg.shuffle_method
        case 'unrestricted'
            % Original v1.0: circular shift across entire session
            r = rand() * session_duration;
            laser_shifted = laser_on + r;
            laser_shifted(laser_shifted > session_duration) = ...
                laser_shifted(laser_shifted > session_duration) - session_duration;

        case 'random_stationary'
            % Randomly place nTrials onset times within stationary periods
            laser_shifted = draw_random_from_intervals(nTrials, stat_iv, stat_durations, stat_duration_total);

        case 'circular_stationary'
            % Circular shift within concatenated stationary timeline
            laser_shifted = circular_shift_in_intervals(laser_on, stat_iv, stat_durations, stat_cumsum, stat_duration_total);
    end

    [disp_shuff, ~] = compute_mean_abs_displacement( ...
        ep_tvec, ep_data_deg, laser_shifted, eye_grid, bl_mask, cfg.window);

    if cfg.smoothEye && cfg.smooth_win > 1
        disp_shuff = smoothdata(disp_shuff, 'gaussian', cfg.smooth_win);
    end

    shuff_matrix(iShuff, :) = disp_shuff;

    if cfg.verbose && mod(iShuff, max(1, floor(cfg.nShuffles / 10))) == 0
        fprintf('    shuffle %d / %d\n', iShuff, cfg.nShuffles);
    end
end

%% ================================================================
%  Statistics
%  ================================================================
shuff_mean = mean(shuff_matrix, 1, 'omitnan');
shuff_sem  = std(shuff_matrix, 0, 1, 'omitnan') / sqrt(cfg.nShuffles);

% Peak |displacement| in test window (one-tailed)
peak_disp_real = max(disp_real_smooth(test_bins));
peak_disp_shuff = NaN(cfg.nShuffles, 1);
for iShuff = 1:cfg.nShuffles
    peak_disp_shuff(iShuff) = max(shuff_matrix(iShuff, test_bins));
end

% Per-bin p-value (one-tailed: fraction of shuffles >= real)
pval_per_bin = NaN(1, nBins);
for iB = 1:nBins
    pval_per_bin(iB) = sum(shuff_matrix(:, iB) >= disp_real_smooth(iB)) / cfg.nShuffles;
end

sig_001 = (pval_per_bin < cfg.alpha_001) & test_bins;
sig_05  = (pval_per_bin < cfg.alpha_05) & ~sig_001 & test_bins;

percentile_overall = sum(peak_disp_shuff < peak_disp_real) / cfg.nShuffles;

if cfg.verbose
    fprintf('  Test window: [%.2f, %.2f] s\n', cfg.test_window(1), cfg.test_window(2));
    fprintf('  Real peak |displacement| in test window: %.3f deg\n', peak_disp_real);
    fprintf('  Shuffle peak |displacement|: %.3f +/- %.3f deg (mean +/- SD)\n', ...
        mean(peak_disp_shuff, 'omitnan'), std(peak_disp_shuff, 'omitnan'));
    fprintf('  Overall percentile: %.4f\n', percentile_overall);
    fprintf('  Bins with p < .001: %d | p < .05: %d\n', sum(sig_001), sum(sig_05));
end

%% Pack output
shuffle_out = struct();
shuffle_out.eye_grid            = eye_grid;
shuffle_out.disp_real           = disp_real;
shuffle_out.disp_real_smooth    = disp_real_smooth;
shuffle_out.shuff_matrix        = shuff_matrix;
shuffle_out.shuff_mean          = shuff_mean;
shuffle_out.shuff_sem           = shuff_sem;
shuffle_out.pval_per_bin        = pval_per_bin;
shuffle_out.sig_001             = sig_001;
shuffle_out.sig_05              = sig_05;
shuffle_out.laser_bins          = laser_bins;
shuffle_out.test_bins           = test_bins;
shuffle_out.test_window         = cfg.test_window;
shuffle_out.peak_disp_real      = peak_disp_real;
shuffle_out.peak_disp_shuff     = peak_disp_shuff;
shuffle_out.percentile_overall  = percentile_overall;
shuffle_out.n_valid_trials      = n_valid_real;
shuffle_out.shuffle_method      = cfg.shuffle_method;
shuffle_out.stat_duration_total = stat_duration_total;

end

%% ================================================================
%  LOCAL: compute mean |displacement| for a set of event times
%  ================================================================
function [mean_abs_disp, n_valid] = compute_mean_abs_displacement( ...
        ep_tvec, ep_data_deg, event_times, eye_grid, bl_mask, window)

    nTrials = length(event_times);
    nBins   = length(eye_grid);
    trial_traces = NaN(nTrials, nBins);

    for iT = 1:nTrials
        t0 = event_times(iT);
        search_mask = ep_tvec >= (t0 + window(1) - 0.1) & ...
                      ep_tvec <= (t0 + window(2) + 0.1);
        if sum(search_mask) < 5; continue; end

        rel_times = ep_tvec(search_mask) - t0;
        rel_data  = ep_data_deg(search_mask);
        interp_trace = interp1(rel_times, rel_data, eye_grid, 'linear', NaN);

        bl_val = mean(interp_trace(bl_mask), 'omitnan');
        if isnan(bl_val); continue; end

        trial_traces(iT, :) = abs(interp_trace - bl_val);
    end

    valid_mask = ~all(isnan(trial_traces), 2);
    n_valid    = sum(valid_mask);

    if n_valid >= 1
        mean_abs_disp = mean(trial_traces(valid_mask, :), 1, 'omitnan');
    else
        mean_abs_disp = NaN(1, nBins);
    end
end

%% ================================================================
%  LOCAL: draw N random times from within stationary intervals
%  ================================================================
function onsets = draw_random_from_intervals(nTrials, stat_iv, stat_durations, stat_total)
% Uniformly sample nTrials times from within the union of stationary intervals.
% Uses inverse CDF: draw uniform [0, stat_total], then map to real time.

    u = rand(nTrials, 1) * stat_total;  % uniform draws in virtual timeline
    onsets = NaN(nTrials, 1);

    cum_dur = cumsum(stat_durations);
    for iT = 1:nTrials
        % Find which interval this draw falls in
        int_idx = find(cum_dur >= u(iT), 1, 'first');
        if int_idx == 1
            offset_in_interval = u(iT);
        else
            offset_in_interval = u(iT) - cum_dur(int_idx - 1);
        end
        onsets(iT) = stat_iv(int_idx, 1) + offset_in_interval;
    end
end

%% ================================================================
%  LOCAL: circular shift within concatenated stationary timeline
%  ================================================================
function laser_shifted = circular_shift_in_intervals(laser_on, stat_iv, stat_durations, stat_cumsum, stat_total)
% Map real laser_on times into a virtual stationary timeline (0 to stat_total),
% apply a single random circular shift, then map back to real time.
% Laser times that fall outside stationary periods are randomly placed.

    nTrials = length(laser_on);

    % Step 1: map each laser_on to virtual time
    virtual_times = NaN(nTrials, 1);
    in_stationary = false(nTrials, 1);
    nIntervals = size(stat_iv, 1);

    for iT = 1:nTrials
        for iInt = 1:nIntervals
            if laser_on(iT) >= stat_iv(iInt, 1) && laser_on(iT) <= stat_iv(iInt, 2)
                offset = laser_on(iT) - stat_iv(iInt, 1);
                if iInt == 1
                    virtual_times(iT) = offset;
                else
                    virtual_times(iT) = stat_cumsum(iInt - 1) + offset;
                end
                in_stationary(iT) = true;
                break
            end
        end
    end

    % For laser pulses outside stationary periods, assign random virtual times
    n_outside = sum(~in_stationary);
    if n_outside > 0
        virtual_times(~in_stationary) = rand(n_outside, 1) * stat_total;
    end

    % Step 2: circular shift in virtual timeline
    shift_amount = rand() * stat_total;
    virtual_shifted = mod(virtual_times + shift_amount, stat_total);

    % Step 3: map back to real time
    laser_shifted = NaN(nTrials, 1);
    for iT = 1:nTrials
        vt = virtual_shifted(iT);
        int_idx = find(stat_cumsum >= vt, 1, 'first');
        if int_idx == 1
            offset_in_interval = vt;
        else
            offset_in_interval = vt - stat_cumsum(int_idx - 1);
        end
        laser_shifted(iT) = stat_iv(int_idx, 1) + offset_in_interval;
    end
end

%% Version History
% 2026-04-07  v1.2  Added cfg.test_window (default [0, 0.2]) — significance
%                    testing now uses a short early window to capture
%                    transient laser-evoked eye movements rather than the
%                    full laser duration. The full trace is still computed
%                    and returned for visualization. Fixes monotonic drift
%                    artifact from abs(baseline-subtracted position).
% 2026-04-07  v1.1  Added stationary-period-restricted shuffling.
%                    Two new methods: 'random_stationary' (randomly place
%                    onsets within stationary intervals) and
%                    'circular_stationary' (circular shift within
%                    concatenated stationary timeline). Default changed
%                    from 'unrestricted' to 'random_stationary'.
%                    Requires cfg.stat_intervals for stationary methods.
% 2026-04-07  v1.0  Initial version. Unrestricted circular shuffle.
