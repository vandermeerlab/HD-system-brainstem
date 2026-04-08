function [results, params] = screen_AHV_vs_EV_opposition(cfg_in)
% screen_AHV_vs_EV_opposition  Screen all NPH neurons for opposite AHV vs slow-phase EV tuning.
%
%   For the AHV vs. EV supplemental figure (Figure S2?), we need a neuron
%   whose AHV tuning curve and slow-phase eye velocity tuning curve are
%   clearly "opposite" — e.g., a neuron that prefers CW head rotation but
%   CCW slow-phase eye velocity. This function loops through all neurons in
%   the master t-file list, computes both TCs, and scores how opposite they are.
%
%   INPUTS:
%       cfg_in              - config structure with optional fields:
%           .tfilelist_path - path to NPH_all_tfilelist.mat
%                             (default: 'C:\Jeff\U01\datatouse\NPH all tfilelist.mat')
%           .AHV_binEdges   - bin edges for AHV TC (default: linspace(-200,200,101))
%           .EV_binEdges    - bin edges for EV TC  (default: linspace(-25,25,100))
%           .AHV_minOcc     - min occupancy for AHV bins (default: 100)
%           .EV_minOcc      - min occupancy for EV bins  (default: 25)
%           .saccade_pre    - seconds to cut before saccade (default: 0.1)
%           .saccade_post   - seconds to cut after saccade  (default: 0.1)
%           .AHV_slope_range - AHV range for slope fit [min max] (default: [-150 150])
%           .EV_slope_range  - EV range for slope fit [min max] (default: [-10 10])
%           .doPlot         - plot top N candidates (default: 1)
%           .nPlot          - how many top candidates to plot (default: 5)
%           .savePath       - where to save results and figs (default: pwd)
%           .doSave         - save results and figures (default: 1)
%
%   OUTPUTS:
%       results             - structure with fields:
%           .tfilelist      - cell array of t-file paths
%           .session_dirs   - cell array of session directories
%           .AHV_slopes     - [nNeurons x 1] slope of AHV TC
%           .EV_slopes      - [nNeurons x 1] slope of slow-phase EV TC
%           .opposition_score - [nNeurons x 1] metric (more negative = more opposite)
%           .AHV_tc         - {nNeurons x 1} cell array of AHV TC structs
%           .EV_tc          - {nNeurons x 1} cell array of EV TC structs
%           .AHV_binCenters - bin centers for AHV
%           .EV_binCenters  - bin centers for EV
%           .skipped        - logical, neurons that were skipped (missing data etc.)
%           .skip_reason    - cell array of skip reasons
%       params              - structure with all parameter values used
%
%   USAGE:
%       cfg = [];
%       cfg.doPlot = 1;
%       cfg.nPlot = 5;
%       [results, params] = screen_AHV_vs_EV_opposition(cfg);
%
%   DEPENDENCIES:
%       LoadSessionData, TuningCurves, SelectTS, antirestrict, HD_GetSSN,
%       MakeQfromS, nearest_idx3, tsd, ProcessConfig2
%       Requires: NPH all tfilelist.mat, *saccades-edited.mat per session
%

%% Config defaults
cfg_def = [];
cfg_def.tfilelist_path = 'C:\Jeff\U01\datatouse\NPH all tfilelist.mat';
cfg_def.AHV_binEdges = {linspace(-200, 200, 101)};
cfg_def.EV_binEdges = {linspace(-40, 40, 81)};   % deg/s; slow-phase EV during VOR typically within +-35
cfg_def.AHV_minOcc = 100;
cfg_def.EV_minOcc = 25;
cfg_def.saccade_pre = 0.1;
cfg_def.saccade_post = 0.1;
cfg_def.AHV_slope_range = [-150 150];
cfg_def.EV_slope_range = [-30 30];   % deg/s; fit slope over central range of slow-phase EV
cfg_def.doPlot = 1;
cfg_def.nPlot = 5;
cfg_def.savePath = 'C:\Jeff\U01\datatouse\___analysis_output';
cfg_def.doSave = 1;
cfg_def.FontSize = 24;
cfg_def.verbose = 1;

cfg = ProcessConfig2(cfg_def, cfg_in);

%% Store params
params = cfg;
params.function_name = 'screen_AHV_vs_EV_opposition';
params.date_run = datestr(now, 'yyyy-mm-dd');

%% Create save folder with date and timestamp
if cfg.doSave
    timestamp_str = datestr(now, 'yyyy-mm-dd_hh-MM-AM');  % e.g., 2026-03-17_02-22-PM
    save_subfolder = fullfile(cfg.savePath, ['screen_AHV_vs_EV_' timestamp_str]);
    if ~exist(save_subfolder, 'dir')
        mkdir(save_subfolder);
    end
    params.save_subfolder = save_subfolder;
end

%% Load master neuron list
load(cfg.tfilelist_path, 'NPH_all_tfilelist');
nNeurons = length(NPH_all_tfilelist);
if cfg.verbose; fprintf('Loaded %d neurons from tfilelist.\n', nNeurons); end

%% Identify unique sessions
% Each t-file path looks like: 'C:\Jeff\U01\...\M222-2022-01-21\M222-2022-01-21-TT01_2.t'
% We group by session directory so we only load session data once per session.
session_dirs_all = cell(nNeurons, 1);
for iN = 1:nNeurons
    session_dirs_all{iN} = fileparts(NPH_all_tfilelist{iN});
end
[unique_sessions, ~, session_idx] = unique(session_dirs_all);
nSessions = length(unique_sessions);
if cfg.verbose; fprintf('Found %d unique sessions.\n', nSessions); end

%% Pre-allocate results
AHV_slopes = NaN(nNeurons, 1);
EV_slopes = NaN(nNeurons, 1);
opposition_scores = NaN(nNeurons, 1);
AHV_tc_all = cell(nNeurons, 1);
EV_tc_all = cell(nNeurons, 1);
skipped = false(nNeurons, 1);
skip_reason = cell(nNeurons, 1);

%% Loop through sessions
for iSess = 1:nSessions
    fd = unique_sessions{iSess};
    neuron_indices = find(session_idx == iSess);  % which neurons belong to this session
    nCellsThisSession = length(neuron_indices);
    
    if cfg.verbose
        fprintf('\n=== Session %d/%d: %s (%d neurons) ===\n', iSess, nSessions, fd, nCellsThisSession);
    end
    
    %% Load session data
    try
        pushdir(fd);
    catch
        for jj = 1:nCellsThisSession
            skipped(neuron_indices(jj)) = true;
            skip_reason{neuron_indices(jj)} = 'cannot pushdir to session folder';
        end
        continue
    end
    
    try
        sd = LoadSessionData([]);
    catch ME
        if cfg.verbose; fprintf('  LoadSessionData failed: %s\n', ME.message); end
        for jj = 1:nCellsThisSession
            skipped(neuron_indices(jj)) = true;
            skip_reason{neuron_indices(jj)} = ['LoadSessionData error: ' ME.message];
        end
        popdir;
        continue
    end
    
    %% Load saccade data for de-saccading
    SSN = HD_GetSSN;
    saccade_file = strcat(SSN, '-saccades-edited.mat');
    if exist(saccade_file, 'file') ~= 2
        if cfg.verbose; fprintf('  Saccade file not found. Skipping session.\n'); end
        for jj = 1:nCellsThisSession
            skipped(neuron_indices(jj)) = true;
            skip_reason{neuron_indices(jj)} = 'saccades-edited.mat not found';
        end
        popdir;
        continue
    end
    
    try
        saccade_data = load(saccade_file, 'combinedSaccades', 'diffHdeg');
        combinedSaccades = saccade_data.combinedSaccades;
        diffHdeg_raw = saccade_data.diffHdeg;
    catch ME
        if cfg.verbose; fprintf('  Error loading saccade file: %s\n', ME.message); end
        for jj = 1:nCellsThisSession
            skipped(neuron_indices(jj)) = true;
            skip_reason{neuron_indices(jj)} = ['saccade file load error: ' ME.message];
        end
        popdir;
        continue
    end
    
    %% De-saccade the AHV and EV signals
    % Get saccade onset times (remove NaNs)
    timestouse = ~isnan(combinedSaccades);
    combinedSaccadesToUse = combinedSaccades(timestouse);
    pre_sacc = combinedSaccadesToUse - cfg.saccade_pre;
    post_sacc = combinedSaccadesToUse + cfg.saccade_post;
    
    % De-saccade AHV
    [~, keepAHV] = restrict(sd.AHV, pre_sacc, post_sacc);
    keepersAHV = keepAHV == 0;  % invert: 1 = keep (outside saccade window)
    AHVr = tsd(sd.AHV.tvec(keepersAHV), sd.AHV.data(keepersAHV));
    
    % De-saccade EV (eye velocity in degrees)
    % Ensure correct orientation for tsd operations
    if size(diffHdeg_raw.tvec, 1) == 1
        diffHdeg_raw.tvec = diffHdeg_raw.tvec';
    end
    if size(diffHdeg_raw.data, 1) ~= 1
        diffHdeg_raw.data = diffHdeg_raw.data';
    end
    
    [~, keepEV] = restrict(diffHdeg_raw, pre_sacc, post_sacc);
    keepersEV = keepEV == 0;
    diffHr = tsd(diffHdeg_raw.tvec(keepersEV), diffHdeg_raw.data(keepersEV));
    
    %% Loop through neurons in this session
    for jj = 1:nCellsThisSession
        iN = neuron_indices(jj);  % index into the master list
        
        % Figure out which cell index this is within the session
        [~, tfile_name] = fileparts(NPH_all_tfilelist{iN});
        % Match to sd.fn
        cell_idx = [];
        for kk = 1:length(sd.fn)
            fn_clean = strrep(sd.fn{kk}, '-', '-');
            tfile_clean = strrep(tfile_name, '_', '-');
            if strcmp(fn_clean, tfile_clean)
                cell_idx = kk;
                break
            end
        end
        
        if isempty(cell_idx)
            skipped(iN) = true;
            skip_reason{iN} = 'could not match t-file to sd.fn';
            continue
        end
        
        try
            %% Select this neuron's spike train and de-saccade it
            myCell = SelectTS([], sd.S, cell_idx);
            [~, keepSpk] = restrict(myCell, pre_sacc, post_sacc);
            myCellr = myCell;
            myCellr.t{1,1}(keepSpk{1,1} == 1) = [];  % remove saccade-window spikes
            
            %% Compute AHV tuning curve (slow-phase only)
            cfg_tcAHV = [];
            cfg_tcAHV.nBins = length(cfg.AHV_binEdges{1}) - 1;
            cfg_tcAHV.binEdges = cfg.AHV_binEdges;
            cfg_tcAHV.minOcc = cfg.AHV_minOcc;
            cfg_tcAHV.occ_dt = median(diff(AHVr.tvec));
            tc_AHV = TuningCurves(cfg_tcAHV, myCellr, AHVr);
            AHV_tc_all{iN} = tc_AHV;
            
            %% Compute EV tuning curve (slow-phase only)
            cfg_tcEV = [];
            cfg_tcEV.nBins = length(cfg.EV_binEdges{1}) - 1;
            cfg_tcEV.binEdges = cfg.EV_binEdges;
            cfg_tcEV.minOcc = cfg.EV_minOcc;
            cfg_tcEV.occ_dt = median(diff(diffHr.tvec));
            tc_EV = TuningCurves(cfg_tcEV, myCellr, diffHr);
            EV_tc_all{iN} = tc_EV;
            
            %% Fit slopes to each TC
            % AHV slope: linear fit over the specified range
            ahv_bc = tc_AHV.binCenters;
            ahv_tc_vals = smoothdata(tc_AHV.tc);
            ahv_range_idx = ahv_bc >= cfg.AHV_slope_range(1) & ahv_bc <= cfg.AHV_slope_range(2);
            ahv_bc_fit = ahv_bc(ahv_range_idx);
            ahv_tc_fit = ahv_tc_vals(ahv_range_idx);
            % Remove NaN bins
            valid_ahv = ~isnan(ahv_tc_fit);
            if sum(valid_ahv) > 5
                p_ahv = polyfit(ahv_bc_fit(valid_ahv), ahv_tc_fit(valid_ahv), 1);
                AHV_slopes(iN) = p_ahv(1);  % slope in Hz per deg/s
            end
            
            % EV slope: linear fit over the specified range
            ev_bc = tc_EV.binCenters;
            ev_tc_vals = smoothdata(tc_EV.tc);
            ev_range_idx = ev_bc >= cfg.EV_slope_range(1) & ev_bc <= cfg.EV_slope_range(2);
            ev_bc_fit = ev_bc(ev_range_idx);
            ev_tc_fit = ev_tc_vals(ev_range_idx);
            valid_ev = ~isnan(ev_tc_fit);
            if sum(valid_ev) > 5
                p_ev = polyfit(ev_bc_fit(valid_ev), ev_tc_fit(valid_ev), 1);
                EV_slopes(iN) = p_ev(1);  % slope in Hz per deg/s
            end
            
            %% Compute opposition score
            % Product of slopes: if AHV slope is negative (CW pref) and EV slope is
            % positive (temporal/CCW pref), the product is negative. More negative = more opposite.
            % We also want both to be steep, so multiply by geometric mean of magnitudes.
            if ~isnan(AHV_slopes(iN)) && ~isnan(EV_slopes(iN))
                opposition_scores(iN) = AHV_slopes(iN) * EV_slopes(iN);
                % More negative = better opposition
            end
            
            if cfg.verbose && mod(jj, 5) == 0
                fprintf('  Neuron %d/%d (master idx %d): AHV slope=%.3f, EV slope=%.3f, score=%.3f\n', ...
                    jj, nCellsThisSession, iN, AHV_slopes(iN), EV_slopes(iN), opposition_scores(iN));
            end
            
        catch ME
            skipped(iN) = true;
            skip_reason{iN} = ['TC computation error: ' ME.message];
            if cfg.verbose; fprintf('  Neuron %d error: %s\n', iN, ME.message); end
        end
    end
    
    popdir;
end

%% Package results
results.tfilelist = NPH_all_tfilelist;
results.session_dirs = session_dirs_all;
results.AHV_slopes = AHV_slopes;
results.EV_slopes = EV_slopes;
results.opposition_score = opposition_scores;
results.AHV_tc = AHV_tc_all;
results.EV_tc = EV_tc_all;
results.AHV_binCenters = tc_AHV.binCenters;  % same for all neurons
results.EV_binCenters = tc_EV.binCenters;
results.skipped = skipped;
results.skip_reason = skip_reason;

%% Report summary
nSkipped = sum(skipped);
nProcessed = nNeurons - nSkipped;
fprintf('\n=== SUMMARY ===\n');
fprintf('Processed: %d/%d neurons\n', nProcessed, nNeurons);
fprintf('Skipped:   %d neurons\n', nSkipped);

% Find best candidates (most negative opposition score = most opposite)
valid_idx = find(~isnan(opposition_scores) & ~skipped);
[sorted_scores, sort_order] = sort(opposition_scores(valid_idx), 'ascend');  % most negative first
best_indices = valid_idx(sort_order);

fprintf('\nTop %d most opposite neurons (AHV slope * EV slope most negative):\n', min(cfg.nPlot, length(best_indices)));
fprintf('%-5s  %-50s  %-10s  %-10s  %-12s\n', 'Rank', 'T-file', 'AHV slope', 'EV slope', 'Score');
fprintf('%s\n', repmat('-', 1, 100));
for ii = 1:min(cfg.nPlot, length(best_indices))
    idx = best_indices(ii);
    [~, fname] = fileparts(NPH_all_tfilelist{idx});
    fprintf('%-5d  %-50s  %-10.4f  %-10.4f  %-12.4f\n', ...
        ii, fname, AHV_slopes(idx), EV_slopes(idx), opposition_scores(idx));
end

results.best_indices = best_indices;
results.sorted_scores = sorted_scores;

%% Plot top candidates
if cfg.doPlot && ~isempty(best_indices)
    nToPlot = min(cfg.nPlot, length(best_indices));
    
    for ii = 1:nToPlot
        idx = best_indices(ii);
        [~, fname] = fileparts(NPH_all_tfilelist{idx});
        
        fig_h = figure('Position', [100 100 1200 500]);
        
        % --- AHV TC ---
        subplot(1, 2, 1); hold on;
        ahv_tc = AHV_tc_all{idx};
        plot(ahv_tc.binCenters, smoothdata(ahv_tc.tc), 'k', 'LineWidth', 4);
        xlabel('AHV (deg/s)', 'FontSize', cfg.FontSize);
        ylabel('FR (Hz)', 'FontSize', cfg.FontSize);
        title(sprintf('AHV TC (slope=%.3f)', AHV_slopes(idx)), 'FontSize', cfg.FontSize);
        set(gca, 'FontSize', cfg.FontSize - 4);
        c_ahv = axis;
        axis([c_ahv(1) c_ahv(2) 0 c_ahv(4)]);
        line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
        % Add CW/CCW labels
        text(0.15, 0.9, 'CW', 'Units', 'normalized', 'FontSize', cfg.FontSize - 2);
        text(0.75, 0.9, 'CCW', 'Units', 'normalized', 'FontSize', cfg.FontSize - 2);
        
        % --- EV TC ---
        subplot(1, 2, 2); hold on;
        ev_tc = EV_tc_all{idx};
        plot(ev_tc.binCenters, smoothdata(ev_tc.tc), 'k', 'LineWidth', 4);
        xlabel('Eye Velocity (deg/s)', 'FontSize', cfg.FontSize);
        ylabel('FR (Hz)', 'FontSize', cfg.FontSize);
        title(sprintf('EV TC (slope=%.3f)', EV_slopes(idx)), 'FontSize', cfg.FontSize);
        set(gca, 'FontSize', cfg.FontSize - 4);
        c_ev = axis;
        axis([c_ev(1) c_ev(2) 0 c_ev(4)]);
        line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
        % Add nasal/temporal labels
        text(0.1, 0.9, 'nasal', 'Units', 'normalized', 'FontSize', cfg.FontSize - 2);
        text(0.7, 0.9, 'temporal', 'Units', 'normalized', 'FontSize', cfg.FontSize - 2);
        
        sgtitle(sprintf('Rank #%d: %s  (score=%.3f)', ii, strrep(fname, '_', '-'), opposition_scores(idx)), ...
            'FontSize', cfg.FontSize, 'Interpreter', 'none');
        
        % Embed metadata in figure
        fig_h.UserData.created_by = 'screen_AHV_vs_EV_opposition';
        fig_h.UserData.date_created = datestr(now, 'yyyy-mm-dd');
        fig_h.UserData.neuron_tfile = NPH_all_tfilelist{idx};
        fig_h.UserData.opposition_score = opposition_scores(idx);
        
        % Save figure
        if cfg.doSave
            fig_filename = sprintf('rank%02d_%s', ii, strrep(fname, '_', '-'));
            savefig(fig_h, fullfile(save_subfolder, [fig_filename '.fig']));
            print(fig_h, fullfile(save_subfolder, [fig_filename '.png']), '-dpng', '-r300');
            if cfg.verbose; fprintf('Saved figure: %s\n', fig_filename); end
        end
    end
end

%% Save results
if cfg.doSave
    save(fullfile(save_subfolder, 'screen_AHV_vs_EV_results.mat'), 'results', 'params', '-v7.3');
    if cfg.verbose; fprintf('\nResults saved to: %s\n', save_subfolder); end
end

end

%% Version History
% 2026-03-17  v1.0 Initial version.
