function [ev_out, params] = compute_slow_phase_EV(varargin)
% compute_slow_phase_EV  Compute slow-phase eye velocity for a single session.
%
%   Standalone EV pipeline using Order B:
%     1. Load raw EV (diffHdeg) and saccade times from saccades-edited.mat
%     2. Load stationary period times from AHV_StationaryTimes.mat
%     3. Desaccade (remove samples within saccade_pre/post window)
%     4. Remove platform-stationary periods
%     5. Threshold on raw desaccaded signal (|EV| < velocity_threshold)
%     6. Median filter (medfilt1)
%     7. Interpolate within contiguous segments onto uniform grid, NaN across gaps
%
%   This replaces the EV-specific portions of getSlowPhaseData2.m with a
%   lighter, self-contained function that uses Order B (desaccade -> threshold
%   -> medfilt1) instead of Order A (medfilt1 -> smoothdata -> threshold).
%
%   INPUTS:
%       varargin    - name-value pairs, e.g.:
%                     compute_slow_phase_EV('velocity_threshold', 40)
%
%   NAME-VALUE OPTIONS:
%       sess_dir            - session directory (default: pwd)
%       saccade_pre         - seconds to exclude before saccade (default: 0.1)
%       saccade_post        - seconds to exclude after saccade (default: 0.1)
%       velocity_threshold  - deg/s slow-phase cutoff (default: 35)
%       medfilt_win         - median filter window in samples (default: 11)
%       gap_thresh          - seconds; gaps larger than this define segment
%                             boundaries for interpolation (default: 0.1)
%       verbose             - print progress messages (default: 1)
%
%   OUTPUTS:
%       ev_out              - structure with fields:
%           .ev_raw         - tsd of raw EV (from diffHdeg, no processing)
%           .ev_gappy       - tsd of slow-phase EV after Order B pipeline
%                             (gappy: samples physically removed)
%           .ev_interp      - tsd on uniform grid with NaN across gaps.
%                             Ready for TuningCurves (NaN bins excluded).
%           .ev_dt          - scalar, time step of uniform grid (seconds)
%           .n_saccades     - number of saccade events used for exclusion
%           .n_stationary   - number of stationary periods excluded
%           .n_samples_raw  - total raw EV samples
%           .n_samples_slow - samples surviving Order B pipeline (gappy)
%           .n_samples_valid- non-NaN samples on interpolated grid
%           .SSN            - session short name
%
%       params              - structure with all parameter values used
%
%   USAGE:
%       % From within session directory:
%       cd('C:\Jeff\U01\datatouse\M222\M222-2022-01-21')
%       [ev_out, params] = compute_slow_phase_EV();
%
%       % Specifying session directory:
%       [ev_out, params] = compute_slow_phase_EV('sess_dir', 'C:\...\M222-2022-01-21');
%
%       % Custom threshold:
%       [ev_out, params] = compute_slow_phase_EV('velocity_threshold', 40, 'medfilt_win', 15);
%
%   DEPENDENCIES:
%       HD_GetSSN, restrict, antirestrict, tsd, process_varargin
%
%                   Extracts Order B EV pipeline from plot_AHV_vs_EV_figure
%                   v1.3 and getSlowPhaseData2. Adds stationary-period
%                   exclusion (via antirestrict) which was missing from
%                   plot_AHV_vs_EV_figure's inline EV computation.
% =========================================================================

%% Defaults
sess_dir            = pwd;
saccade_pre         = 0.1;       % s before saccade to exclude
saccade_post        = 0.1;       % s after saccade to exclude
velocity_threshold  = 35;        % deg/s
medfilt_win         = 11;        % samples
gap_thresh          = 0.1;       % s — segment boundary threshold for interpolation
remove_stationary   = 0;         % 0 = keep stationary periods, 1 = remove them
verbose             = 1;

%% Process varargin
process_varargin(varargin);

%% Store params
params = struct();
params.function_name        = 'compute_slow_phase_EV';
params.version              = '1.0';
params.date_run             = datestr(now, 'yyyy-mm-dd');
params.sess_dir             = sess_dir;
params.saccade_pre          = saccade_pre;
params.saccade_post         = saccade_post;
params.velocity_threshold   = velocity_threshold;
params.medfilt_win          = medfilt_win;
params.gap_thresh           = gap_thresh;
params.remove_stationary    = remove_stationary;
params.pipeline             = 'Order B: desaccade -> [optional: remove stationary] -> threshold -> medfilt1 -> interpolate';

%% Navigate to session
pushdir(sess_dir);
SSN = HD_GetSSN;
if verbose; fprintf('compute_slow_phase_EV: %s\n', SSN); end

%% Load data
% Saccade file: raw EV and curated saccade times
sacc_file = strcat(SSN, '-saccades-edited.mat');
if ~exist(sacc_file, 'file')
    error('Saccade file not found: %s', sacc_file);
end
sacc_data = load(sacc_file, 'combinedSaccades', 'diffHdeg');
diffHdeg_raw    = sacc_data.diffHdeg;
combinedSacc    = sacc_data.combinedSaccades;

% Stationary times file (only needed if removing stationary periods)
if remove_stationary
    stat_file = strcat(SSN, '-AHV_StationaryTimes.mat');
    if ~exist(stat_file, 'file')
        error('Stationary times file not found: %s', stat_file);
    end
    stat_data = load(stat_file, 'STtstart', 'STtend');
    STtstart = stat_data.STtstart;
    STtend   = stat_data.STtend;
end

%% Ensure correct orientation (column tvec, row data)
if size(diffHdeg_raw.tvec, 1) == 1
    diffHdeg_raw.tvec = diffHdeg_raw.tvec';
end
if size(diffHdeg_raw.data, 1) ~= 1
    diffHdeg_raw.data = diffHdeg_raw.data';
end

n_samples_raw = length(diffHdeg_raw.tvec);
if verbose; fprintf('  Raw EV: %d samples\n', n_samples_raw); end

%% Store raw EV (before any processing)
ev_raw_tsd = tsd(diffHdeg_raw.tvec, diffHdeg_raw.data);

%% Step 1: Desaccade — remove samples within saccade windows
sacc_times  = combinedSacc(~isnan(combinedSacc));
pre_sacc    = sacc_times - saccade_pre;
post_sacc   = sacc_times + saccade_post;
n_saccades  = length(sacc_times);

[~, keep_mask_sacc] = restrict(diffHdeg_raw, pre_sacc, post_sacc);
is_outside_saccade  = keep_mask_sacc == 0;  % invert: keep samples OUTSIDE saccade windows

tvec_desacc = diffHdeg_raw.tvec(is_outside_saccade);
ev_desacc   = diffHdeg_raw.data(is_outside_saccade);

if verbose; fprintf('  After desaccade: %d samples (%d saccades excluded)\n', length(tvec_desacc), n_saccades); end

%% Step 2: Optionally remove stationary periods
if remove_stationary
    desacc_tsd = tsd(tvec_desacc, ev_desacc);
    ev_moving  = antirestrict(desacc_tsd, STtstart, STtend);

    tvec_moving = ev_moving.tvec;
    ev_mov_data = ev_moving.data;
    n_stationary = length(STtstart);

    if verbose; fprintf('  After removing stationary: %d samples (%d stationary periods)\n', length(tvec_moving), n_stationary); end
else
    tvec_moving = tvec_desacc;
    ev_mov_data = ev_desacc;
    n_stationary = 0;
    if verbose; fprintf('  Stationary periods: kept (remove_stationary = 0)\n'); end
end

%% Step 3: Threshold on raw desaccaded+moving signal
slow_idx    = abs(ev_mov_data) < velocity_threshold;
tvec_slow   = tvec_moving(slow_idx);
ev_slow     = ev_mov_data(slow_idx);

if verbose; fprintf('  After velocity threshold (|EV| < %d deg/s): %d samples\n', velocity_threshold, length(tvec_slow)); end

%% Step 4: Median filter
ev_filt = medfilt1(ev_slow, medfilt_win);

%% Package gappy tsd (samples physically removed, not on uniform grid)
ev_gappy = tsd(tvec_slow, ev_filt);

n_samples_slow = length(tvec_slow);

%% Step 5: Interpolate within contiguous segments onto uniform grid
ev_dt           = median(diff(tvec_slow));
tvec_gaps       = diff(tvec_slow);
seg_breaks      = find(tvec_gaps > gap_thresh);
seg_starts      = [1; seg_breaks + 1];
seg_ends        = [seg_breaks; length(tvec_slow)];

tvec_uniform    = (tvec_slow(1) : ev_dt : tvec_slow(end))';
ev_uniform      = nan(size(tvec_uniform));

for iSeg = 1:length(seg_starts)
    idx_s   = seg_starts(iSeg);
    idx_e   = seg_ends(iSeg);
    if idx_e - idx_s < 2
        continue  % need at least 3 samples for meaningful interpolation
    end
    t_s     = tvec_slow(idx_s);
    t_e     = tvec_slow(idx_e);
    in_seg  = tvec_uniform >= t_s & tvec_uniform <= t_e;
    ev_uniform(in_seg) = interp1(tvec_slow(idx_s:idx_e), ...
                                  ev_filt(idx_s:idx_e), ...
                                  tvec_uniform(in_seg), 'linear');
end

ev_interp = tsd(tvec_uniform, ev_uniform');

n_samples_valid = sum(~isnan(ev_uniform));

if verbose
    fprintf('  Uniform grid: %d total samples, %d valid (non-NaN), %d segments\n', ...
        length(tvec_uniform), n_samples_valid, length(seg_starts));
    fprintf('  dt = %.4f s (%.1f Hz)\n', ev_dt, 1/ev_dt);
end

%% Package output
ev_out = struct();
ev_out.ev_raw           = ev_raw_tsd;
ev_out.ev_gappy         = ev_gappy;
ev_out.ev_interp        = ev_interp;
ev_out.ev_dt            = ev_dt;
ev_out.n_saccades       = n_saccades;
ev_out.n_stationary     = n_stationary;
ev_out.n_samples_raw    = n_samples_raw;
ev_out.n_samples_slow   = n_samples_slow;
ev_out.n_samples_valid  = n_samples_valid;
ev_out.SSN              = SSN;

popdir;

if verbose; fprintf('  Done.\n'); end

end

%% Version History
% 2026-03-17  v1.0 Initial version.
