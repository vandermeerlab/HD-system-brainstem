function [shuffle_out, params] = compute_optotag_shuffle(spike_times, laser_on, session_duration, cfg_in, varargin)
% compute_optotag_shuffle  Circular shuffle significance test for optotag PETHs
%
%   Circularly shifts laser onset times relative to a fixed spike train to
%   generate a null distribution of PETHs. Tests whether real laser-evoked
%   firing rate exceeds chance at two levels:
%     (1) Overall: mean FR in laser window vs. shuffle distribution
%     (2) Per-bin: FR in each time bin vs. shuffle distribution at that bin
%
%   INPUTS:
%       spike_times      - [nSpikes x 1] spike timestamps (session-relative, s)
%       laser_on         - [nTrials x 1] laser onset timestamps (session-relative, s)
%       session_duration - scalar, total session duration in seconds
%       cfg_in           - config structure (all fields optional)
%       varargin         - name-value pairs overriding any cfg field
%
%   CFG FIELDS:
%       .nShuffles      - number of circular shuffles (default: 10)
%       .window         - PETH window in s (default: [-1 2])
%       .dt             - PETH bin size in s (default: 0.05)
%       .laser_dur      - laser pulse duration in s (default: 1.0)
%       .smoothFR       - smooth FR traces (default: 1)
%       .smooth_win     - Gaussian smoothing window in bins (default: 3)
%       .alpha_001      - threshold for p < .001 (default: 0.001)
%       .alpha_05       - threshold for p < .05 (default: 0.05)
%       .chenhao_method - if 1, pool all shuffle bins into one flat null
%                         distribution (default: 0)
%       .verbose        - print progress (default: 1)
%
%   OUTPUTS:
%       shuffle_out     - struct with fields:
%           .bin_centers       - [1 x nBins] PETH bin centers (s)
%           .fr_real           - [1 x nBins] real mean FR PETH (Hz)
%           .fr_real_smooth    - [1 x nBins] smoothed real FR PETH (Hz)
%           .shuff_matrix      - [nShuffles x nBins] all shuffled FR PETHs
%           .shuff_mean        - [1 x nBins] mean of shuffled PETHs per bin
%           .shuff_sem         - [1 x nBins] SEM of shuffled PETHs per bin
%           .pval_per_bin      - [1 x nBins] per-bin p-value (fraction of
%                                shuffles with FR >= real FR)
%           .sig_001           - [1 x nBins] logical, bins with p < .001
%           .sig_05            - [1 x nBins] logical, bins with p < .05
%                                (but not p < .001)
%           .laser_bins        - [1 x nBins] logical, bins within laser window
%           .peak_FR_real      - scalar, max FR in laser window (real PETH)
%           .peak_FR_shuff     - [nShuffles x 1] max FR in laser window
%                                per shuffle
%           .percentile_overall - scalar, fraction of shuffles with peak FR
%                                < real peak FR (>= 0.95 means sig at p<.05)
%           .shift_values      - [nShuffles x 1] random shift amounts used
%       params          - parameter struct
%
%   DEPENDENCIES:
%       ProcessConfig2, process_varargin
%

%% Config defaults
cfg_def = [];
cfg_def.nShuffles   = 10;
cfg_def.window      = [-1 2];
cfg_def.dt          = 0.05;
cfg_def.laser_dur   = 1.0;
cfg_def.smoothFR    = 1;
cfg_def.smooth_win  = 3;
cfg_def.alpha_001   = 0.001;
cfg_def.alpha_05    = 0.05;
cfg_def.chenhao_method = 0;  % if 1, pool all shuffle bins into one flat null distribution
cfg_def.verbose     = 1;

if nargin < 4 || isempty(cfg_in); cfg_in = []; end
cfg = ProcessConfig2(cfg_def, cfg_in);
process_varargin(varargin);

%% Store params
params                  = cfg;
params.function_name    = 'compute_optotag_shuffle';
params.version          = '1.5';
params.date_run         = datestr(now, 'yyyy-mm-dd');
params.nSpikes          = length(spike_times);
params.nTrials          = length(laser_on);
params.session_duration = session_duration;

%% Set up PETH bins
bin_edges   = cfg.window(1) : cfg.dt : cfg.window(2);
bin_centers = bin_edges(1:end-1) + cfg.dt / 2;
nBins       = length(bin_centers);
nTrials     = length(laser_on);

%% Identify which bins fall within the laser-on window
laser_bins = (bin_centers >= 0) & (bin_centers < cfg.laser_dur);

%% Compute the REAL PETH
fr_hist_real = zeros(1, nBins);
for iT = 1:nTrials
    t0 = laser_on(iT);
    spk_in_win = spike_times(spike_times >= t0 + cfg.window(1) & ...
                              spike_times <  t0 + cfg.window(2));
    spk_rel = spk_in_win - t0;
    if ~isempty(spk_rel)
        fr_hist_real = fr_hist_real + histcounts(spk_rel, bin_edges);
    end
end
fr_real = fr_hist_real / (nTrials * cfg.dt);

if cfg.smoothFR && cfg.smooth_win > 1
    fr_real_smooth = smoothdata(fr_real, 'gaussian', cfg.smooth_win);
else
    fr_real_smooth = fr_real;
end

%% Circular shuffle loop
shuff_matrix        = NaN(cfg.nShuffles, nBins);
shift_values        = NaN(cfg.nShuffles, 1);

if cfg.verbose
    fprintf('compute_optotag_shuffle: running %d shuffles...\n', cfg.nShuffles);
end

for iShuff = 1:cfg.nShuffles
    % Pick a random shift amount uniformly from (0, session_duration)
    r = rand() * session_duration;
    shift_values(iShuff) = r;

    % Shift laser onset times
    laser_shifted = laser_on + r;

    % Circular wrap: any shifted time > session_duration wraps around
    laser_shifted(laser_shifted > session_duration) = ...
        laser_shifted(laser_shifted > session_duration) - session_duration;

    % Compute shuffled PETH
    fr_hist_shuff = zeros(1, nBins);
    for iT = 1:nTrials
        t0 = laser_shifted(iT);
        spk_in_win = spike_times(spike_times >= t0 + cfg.window(1) & ...
                                  spike_times <  t0 + cfg.window(2));
        spk_rel_shuff = spk_in_win - t0;
        if ~isempty(spk_rel_shuff)
            fr_hist_shuff = fr_hist_shuff + histcounts(spk_rel_shuff, bin_edges);
        end
    end
    shuff_matrix(iShuff, :) = fr_hist_shuff / (nTrials * cfg.dt);

    if cfg.verbose && mod(iShuff, max(1, floor(cfg.nShuffles / 10))) == 0
        fprintf('  shuffle %d / %d\n', iShuff, cfg.nShuffles);
    end
end

%% Apply identical smoothing to shuffle matrix
if cfg.smoothFR && cfg.smooth_win > 1
    for iShuff = 1:cfg.nShuffles
        shuff_matrix(iShuff, :) = smoothdata(shuff_matrix(iShuff, :), 'gaussian', cfg.smooth_win);
    end
end

%% Shuffle statistics: mean and SEM per bin (on smoothed data)
if cfg.chenhao_method
    % Chenhao method: pool ALL shuffle bins into one flat null distribution.
    % Under the null, every time bin is equivalent, so use a single global
    % mean and SEM. Gives a flat line + constant band for plotting, and
    % per-bin significance is tested against the pooled distribution.
    all_shuff_vals  = shuff_matrix(:);  % pool all nShuffles x nBins values
    pooled_mean     = mean(all_shuff_vals);
    pooled_std      = std(all_shuff_vals);
    pooled_sem      = pooled_std / sqrt(numel(all_shuff_vals));
    shuff_mean      = repmat(pooled_mean, 1, nBins);
    shuff_sem       = repmat(pooled_sem, 1, nBins);
else
    % Per-bin method: independent mean and SEM at each time bin
    shuff_mean = mean(shuff_matrix, 1);
    shuff_sem  = std(shuff_matrix, 0, 1) / sqrt(cfg.nShuffles);
end

%% Peak FR in laser window (one-tailed test, using smoothed traces)
% Real: max smoothed FR bin within laser window
% Each shuffle: max smoothed FR bin within laser window
% Tests whether laser drove spiking above what random alignment produces.
peak_FR_real = max(fr_real_smooth(laser_bins));
peak_FR_shuff = NaN(cfg.nShuffles, 1);
for iShuff = 1:cfg.nShuffles
    peak_FR_shuff(iShuff) = max(shuff_matrix(iShuff, laser_bins));
end

%% Per-bin significance (one-tailed, on smoothed traces)
pval_per_bin = NaN(1, nBins);
if cfg.chenhao_method
    % Test each bin against the pooled null distribution
    for iB = 1:nBins
        pval_per_bin(iB) = sum(all_shuff_vals >= fr_real_smooth(iB)) / numel(all_shuff_vals);
    end
else
    % Test each bin against its own shuffle distribution
    for iB = 1:nBins
        pval_per_bin(iB) = sum(shuff_matrix(:, iB) >= fr_real_smooth(iB)) / cfg.nShuffles;
    end
end

% Significance masks (only within laser window)
sig_001 = (pval_per_bin < cfg.alpha_001) & laser_bins;
sig_05  = (pval_per_bin < cfg.alpha_05) & ~sig_001 & laser_bins;

%% Overall significance (peak FR in laser window, one-tailed)
% Percentile rank: fraction of shuffles with peak FR < real peak FR
% High percentile (e.g., >= 0.95) means significant
percentile_overall = sum(peak_FR_shuff < peak_FR_real) / cfg.nShuffles;

if cfg.verbose
    if cfg.chenhao_method
        fprintf('  Method: chenhao (pooled null distribution)\n');
    else
        fprintf('  Method: per-bin null distribution\n');
    end
    fprintf('  Real peak FR in laser window: %.2f Hz\n', peak_FR_real);
    fprintf('  Shuffle peak FR: %.2f +/- %.2f Hz (mean +/- SD)\n', ...
        mean(peak_FR_shuff), std(peak_FR_shuff));
    fprintf('  Overall percentile: %.4f\n', percentile_overall);
    fprintf('  Bins with p < .001: %d\n', sum(sig_001));
    fprintf('  Bins with p < .05 (not .001): %d\n', sum(sig_05));
end

%% Pack output
shuffle_out = struct();
shuffle_out.bin_centers     = bin_centers;
shuffle_out.fr_real         = fr_real;
shuffle_out.fr_real_smooth  = fr_real_smooth;
shuffle_out.shuff_matrix    = shuff_matrix;
shuffle_out.shuff_mean      = shuff_mean;
shuffle_out.shuff_sem       = shuff_sem;
shuffle_out.pval_per_bin    = pval_per_bin;
shuffle_out.sig_001         = sig_001;
shuffle_out.sig_05          = sig_05;
shuffle_out.laser_bins          = laser_bins;
shuffle_out.peak_FR_real        = peak_FR_real;
shuffle_out.peak_FR_shuff       = peak_FR_shuff;
shuffle_out.percentile_overall  = percentile_overall;
shuffle_out.shift_values        = shift_values;

end

%% Version History
% 2026-03-18  v1.5  Added cfg.chenhao_method option: pools all shuffle bins
%                    into one flat null distribution for per-bin testing and
%                    plotting. Gives flat shuffle mean + constant SEM band.
% 2026-03-18  v1.4  Apply identical Gaussian smoothing to both real and
%                    shuffled PETHs before all significance tests (per-bin
%                    and overall). Eliminates fragmented significance bars
%                    caused by bin-to-bin noise in unsmoothed data.
% 2026-03-18  v1.3  Switched overall test statistic to peak FR in laser
%                    window (one-tailed). Tests whether laser drove spiking
%                    above chance. Simpler and more appropriate for ChR2
%                    excitation than two-tailed deviation tests.
% 2026-03-18  v1.2  Fixed overall test statistic: each PETH's max deviation
%                    is now computed against its OWN mean (matching
%                    saccade_sig_distr_3 logic), not the grand shuffle mean.
%                    Switched output to percentile convention (high = sig).
% 2026-03-18  v1.1  Switched overall test statistic from mean FR in laser
%                    window to max absolute deviation from shuffle mean
%                    (matches saccade_sig_distr_3 approach). Catches
%                    transient responses that cancel out in a window mean.
% 2026-03-18  v1.0  Initial version. Circular shuffle of laser onset times,
%                    per-bin and overall significance tests.
