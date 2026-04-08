function [fig_handles, params, neuron_data] = plot_AHV_vs_EV_figure(cfg_in, varargin)
% plot_AHV_vs_EV_figure  Generate Figure S2: AHV and EV are inversely related.
%
%   Produces a 3-panel figure for each specified neuron:
%     Panel A: Behavioral scatterplot of AHV vs slow-phase EV (session-level)
%     Panel B: AHV tuning curve with raw firing rate data in grey
%     Panel C: Slow-phase EV tuning curve with raw firing rate data in grey
%
%   Panels B and C share the same FR y-axis range.
%
%   EV pipeline via compute_slow_phase_EV (Order B):
%     desaccade -> remove stationary -> threshold -> medfilt1 -> interpolate
%
%   AHV is also desaccaded and has stationary periods removed.
%
%   INPUTS:
%       cfg_in              - config structure (all fields optional)
%       varargin            - name-value pairs overriding any cfg field
%
%   CFG FIELDS:
%       .neuron_indices     - indices into NPH_all_tfilelist to plot
%                             (default: top N from screening results)
%       .neuron_data        - cell array of pre-computed neuron data structs
%                             (from a previous run's 3rd output). If provided,
%                             skips all loading/computation and goes straight
%                             to plotting. Overrides cache_path and results_path.
%       .nPlot              - how many top candidates to plot (default: 5)
%       .results_path       - path to screen_AHV_vs_EV_results.mat
%       .cache_path         - path to previously saved figS2_neuron_data.mat
%       .tfilelist_path     - path to NPH_all_tfilelist.mat
%       .AHV_binEdges       - bin edges for AHV TC (default: linspace(-175,175,89))
%       .EV_binEdges        - bin edges for EV TC  (default: linspace(-40,40,81))
%       .AHV_minOcc         - min occupancy in seconds for AHV bins (default: 1)
%       .EV_minOcc          - min occupancy in seconds for EV bins  (default: 0.5)
%       .saccade_pre        - seconds to cut before saccade (default: 0.1)
%       .saccade_post       - seconds to cut after saccade  (default: 0.1)
%       .doShow             - display figures (default: 0)
%       .doSavePNG          - save .png files (default: 1)
%       .doSaveFIG          - save .fig files (default: 0)
%       .savePath           - root save path
%       .FontSize           - base font size (default: 24)
%       .MarkerSize         - grey raw data dot size (default: 5)
%       .LineWidth          - TC overlay line width (default: 4)
%       .smoothTC           - smooth tuning curves (default: 1)
%       .figWidth           - figure width in pixels (default: 1800)
%       .figHeight          - figure height in pixels (default: 600)
%       .Q_gausswin_sd      - SD for Gaussian FR smoothing in s (default: 0.05)
%       .robustFit          - 0 = regress (OLS), 1 = fitlm with RobustOpts
%                             (bisquare weights, downweights outliers) (default: 0)
%       .panelA_AHV_lim     - Panel A x-axis limits (default: [-175 175])
%       .panelA_EV_lim      - Panel A y-axis limits (default: [-50 50])
%       .panelC_EV_lim      - Panel C x-axis limits (default: [-40 40])
%       .verbose            - print progress (default: 1)
%
%   OUTPUTS:
%       fig_handles         - array of figure handles (one per neuron)
%       params              - structure with all parameter values used
%       neuron_data         - cell array of per-neuron data structs (pass
%                             back in via cfg.neuron_data to skip recompute)
%
%   USAGE:
%       % First run (computes from scratch):
%       cfg = [];
%       cfg.results_path = 'C:\...\screen_AHV_vs_EV_results.mat';
%       [fh, params, nd] = plot_AHV_vs_EV_figure(cfg);
%
%       % Subsequent runs (workspace mode, no recompute):
%       cfg2 = [];
%       cfg2.neuron_data = nd;
%       cfg2.doShow = 1;
%       [fh2, params2] = plot_AHV_vs_EV_figure(cfg2);
%
%   DEPENDENCIES:
%       compute_slow_phase_EV, LoadSessionData, TuningCurves, SelectTS,
%       HD_GetSSN, MakeQfromS, nearest_idx3, tsd, ProcessConfig2,
%       restrict, antirestrict, pushdir, popdir, process_varargin
%

%% Config defaults
cfg_def = [];
cfg_def.neuron_indices          = [];
cfg_def.neuron_data             = {};       % pass in pre-computed data to skip recompute
cfg_def.nPlot                   = 5;
cfg_def.results_path            = '';
cfg_def.cache_path              = '';
cfg_def.tfilelist_path          = 'C:\Jeff\U01\datatouse\NPH all tfilelist.mat';
cfg_def.AHV_binEdges            = {linspace(-175, 175, 89)};  % match panelA_AHV_lim
cfg_def.EV_binEdges             = {linspace(-40, 40, 81)};
cfg_def.AHV_minOcc              = 1;        % seconds (lowered to recover edge bins)
cfg_def.EV_minOcc               = 0.5;      % seconds
cfg_def.saccade_pre             = 0.1;
cfg_def.saccade_post            = 0.1;
cfg_def.doShow                  = 0;
cfg_def.doSavePNG               = 1;
cfg_def.doSaveFIG               = 0;
cfg_def.doSaveEPS               = 0;
cfg_def.savePath                = 'C:\Jeff\U01\datatouse\___analysis_output';
cfg_def.save_subfolder          = '';        % if provided, use this folder (skip auto-generate)
cfg_def.FontSize                = 24;
cfg_def.MarkerSize              = 5;
cfg_def.LineWidth               = 4;
cfg_def.smoothTC                = 1;
cfg_def.smoothTC_win            = 11;       % smoothdata Gaussian window (bins)
cfg_def.figWidth                = 2100;
cfg_def.figHeight               = 650;
cfg_def.panelWidthScale         = 1.15;     % multiplier on default subplot width (1.0 = default)
cfg_def.greyColor               = [0.5 0.5 0.5];
cfg_def.panelA_AHV_lim          = [-175 175];
cfg_def.panelA_EV_lim           = [-50 50];
cfg_def.panelC_EV_lim           = [-40 40];
cfg_def.panelA_useGappy         = 0;        % 0 = interpolated EV (default), 1 = gappy (non-interpolated)
cfg_def.Q_gausswin_sd           = 0.05;
cfg_def.robustFit               = 0;        % 0 = regress (OLS), 1 = fitlm with RobustOpts
cfg_def.remove_stationary       = 0;        % 0 = keep stationary periods, 1 = remove via antirestrict
cfg_def.verbose                 = 1;

cfg = ProcessConfig2(cfg_def, cfg_in);
process_varargin(varargin);  % allow name-value overrides

%% Store params
params                  = cfg;
params.function_name    = 'plot_AHV_vs_EV_figure';
params.version          = '2.3';
params.date_run         = datestr(now, 'yyyy-mm-dd');
params.EV_pipeline      = 'compute_slow_phase_EV (Order B with stationary removal)';

%% Create or use save folder
if isfield(cfg, 'save_subfolder') && ~isempty(cfg.save_subfolder)
    save_subfolder = cfg.save_subfolder;
else
    timestamp_str  = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    save_subfolder = fullfile(cfg.savePath, ['figS2_AHV_vs_EV_' timestamp_str]);
end
if (cfg.doSavePNG || cfg.doSaveFIG || cfg.doSaveEPS) && ~exist(save_subfolder, 'dir')
    mkdir(save_subfolder);
end
params.save_subfolder = save_subfolder;

%% Load from workspace, cache, or compute
if ~isempty(cfg.neuron_data)
    %% Workspace mode — use pre-computed neuron_data passed in via cfg
    neuron_data = cfg.neuron_data;
    nToPlot     = length(neuron_data);
    plot_indices = 1:nToPlot;  % placeholder — indices not meaningful in workspace mode

    % Re-run regression if robustFit setting differs from stored data
    for jj = 1:nToPlot
        nd_j = neuron_data{jj};
        if isempty(nd_j); continue; end

        stored_is_robust = strcmp(nd_j.reg_method, 'fitlm_robust');
        if cfg.robustFit ~= stored_is_robust
            % Recompute regression from stored scatter data
            valid_jj    = ~isnan(nd_j.AHV_for_scatter) & ~isnan(nd_j.EV_for_scatter);
            ahv_col     = nd_j.AHV_for_scatter(valid_jj)';
            ev_col      = nd_j.EV_for_scatter(valid_jj)';

            if cfg.robustFit
                lm_j = fitlm(ahv_col, ev_col, 'RobustOpts', 'on');
                nd_j.reg_slope  = lm_j.Coefficients.Estimate(2);
                nd_j.reg_intcpt = lm_j.Coefficients.Estimate(1);
                nd_j.reg_method = 'fitlm_robust';
                nd_j.reg_R2     = lm_j.Rsquared.Ordinary;
            else
                X_des   = [ones(length(ahv_col), 1), ahv_col];
                [b_j, ~, ~, ~, stats_j] = regress(ev_col, X_des);
                nd_j.reg_slope  = b_j(2);
                nd_j.reg_intcpt = b_j(1);
                nd_j.reg_method = 'regress_OLS';
                nd_j.reg_R2     = stats_j(1);
            end
            neuron_data{jj} = nd_j;

            if cfg.verbose
                fprintf('  Re-ran regression for neuron %d -> %s: slope=%.4f  intercept=%.4f  R2=%.4f\n', ...
                    jj, nd_j.reg_method, nd_j.reg_slope, nd_j.reg_intcpt, nd_j.reg_R2);
            end
        end
    end

    if cfg.verbose
        fprintf('Using pre-computed neuron_data from workspace (%d neurons).\n', nToPlot);
    end

elseif ~isempty(cfg.cache_path) && exist(cfg.cache_path, 'file')
    cached          = load(cfg.cache_path, 'neuron_data', 'plot_indices');
    neuron_data     = cached.neuron_data;
    plot_indices    = cached.plot_indices;
    nToPlot         = length(plot_indices);
    if cfg.verbose
        fprintf('Loaded cached data for %d neurons from:\n  %s\n', nToPlot, cfg.cache_path);
    end

else
    %% Determine which neurons to plot
    load(cfg.tfilelist_path, 'NPH_all_tfilelist');

    if ~isempty(cfg.neuron_indices)
        plot_indices = cfg.neuron_indices;
    elseif ~isempty(cfg.results_path)
        screening       = load(cfg.results_path, 'results');
        plot_indices    = screening.results.best_indices(1:min(cfg.nPlot, ...
                            length(screening.results.best_indices)));
    else
        error('Must provide cfg.neuron_indices, cfg.results_path, or cfg.cache_path.');
    end

    nToPlot     = length(plot_indices);
    neuron_data = cell(nToPlot, 1);

    if cfg.verbose
        fprintf('Computing data for %d neurons.\n', nToPlot);
    end

    %% Loop: load and compute per neuron
    for ii = 1:nToPlot
        iN              = plot_indices(ii);
        tfile_path      = NPH_all_tfilelist{iN};
        [sess_dir, tfile_name] = fileparts(tfile_path);

        if cfg.verbose
            fprintf('\n=== Neuron %d/%d: %s ===\n', ii, nToPlot, tfile_name);
        end

        pushdir(sess_dir);
        sd  = LoadSessionData([]);
        SSN = HD_GetSSN;

        %% Match t-file to sd.fn index
        cell_idx    = [];
        tfile_clean = strrep(tfile_name, '_', '-');
        for kk = 1:length(sd.fn)
            if strcmp(sd.fn{kk}, tfile_clean)
                cell_idx = kk;
                break
            end
        end
        if isempty(cell_idx)
            warning('Could not match %s to sd.fn. Skipping.', tfile_name);
            popdir;
            continue
        end

        %% Load saccade data
        sacc_data           = load(strcat(SSN, '-saccades-edited.mat'), ...
                                'combinedSaccades');
        combinedSaccades    = sacc_data.combinedSaccades;

        %% Build saccade exclusion windows
        sacc_times  = combinedSaccades(~isnan(combinedSaccades));
        pre_sacc    = sacc_times - cfg.saccade_pre;
        post_sacc   = sacc_times + cfg.saccade_post;

        %% Load stationary period times (if needed)
        if cfg.remove_stationary
            stat_file = strcat(SSN, '-AHV_StationaryTimes.mat');
            stat_data = load(stat_file, 'STtstart', 'STtend');
        end

        %% ============================================================
        %  STEP 1: Compute FR from RAW spike train on full AHV timebase
        %  ============================================================
        myCell      = SelectTS([], sd.S, cell_idx);
        occ_dt_AHV  = median(diff(sd.AHV.tvec));

        cfg_Q_ahv               = [];
        cfg_Q_ahv.smooth        = 'gauss';
        cfg_Q_ahv.gausswin_sd   = cfg.Q_gausswin_sd;
        cfg_Q_ahv.dt            = occ_dt_AHV;
        cfg_Q_ahv.tvec_edges    = sd.AHV.tvec(1):occ_dt_AHV:sd.AHV.tvec(end);
        F_ahv_full              = MakeQfromS(cfg_Q_ahv, myCell);
        F_ahv_full.data         = F_ahv_full.data ./ cfg_Q_ahv.dt;

        % Align FR to full AHV timebase
        F_idx_full  = nearest_idx3(sd.AHV.tvec, F_ahv_full.tvec);
        FR_ahv_full = F_ahv_full.data(1, F_idx_full);  % 1×N, same length as sd.AHV.tvec

        %% ============================================================
        %  STEP 2: Restrict AHV + FR together (desaccade, then stationary)
        %  ============================================================
        % Desaccade
        [~, keepAHV]    = restrict(sd.AHV, pre_sacc, post_sacc);
        keepersAHV      = keepAHV == 0;
        AHVr            = tsd(sd.AHV.tvec(keepersAHV), sd.AHV.data(keepersAHV));
        AHV_FR          = FR_ahv_full(keepersAHV);

        if cfg.verbose
            fprintf('  AHV: %d raw -> %d desaccaded\n', ...
                length(sd.AHV.tvec), sum(keepersAHV));
        end

        % Remove stationary periods (if requested)
        if cfg.remove_stationary
            [~, statKeep] = restrict(AHVr, stat_data.STtstart, stat_data.STtend);
            statKeepers   = statKeep == 0;
            AHV_FR        = AHV_FR(statKeepers);
            AHVr          = tsd(AHVr.tvec(statKeepers), AHVr.data(statKeepers));
            if cfg.verbose
                fprintf('  AHV after removing stationary: %d samples\n', length(AHVr.tvec));
            end
        end

        %% ============================================================
        %  STEP 3: Compute slow-phase EV (handles its own desaccade etc.)
        %  ============================================================
        [ev_out, ev_params] = compute_slow_phase_EV('sess_dir', sess_dir, ...
            'saccade_pre', cfg.saccade_pre, 'saccade_post', cfg.saccade_post, ...
            'remove_stationary', cfg.remove_stationary, ...
            'verbose', 0);

        tvec_ev_cont    = ev_out.ev_interp.tvec;
        ev_cont_data    = ev_out.ev_interp.data;
        ev_dt           = ev_out.ev_dt;

        % Ensure correct orientation
        if size(tvec_ev_cont, 1) == 1; tvec_ev_cont = tvec_ev_cont'; end
        if size(ev_cont_data, 1) ~= 1; ev_cont_data = ev_cont_data'; end

        % tsd on uniform grid (for TuningCurves)
        diffHr_slow = tsd(tvec_ev_cont, ev_cont_data);

        %% ============================================================
        %  STEP 4: Compute FR from RAW spike train on EV timebase
        %  ============================================================
        cfg_Q_ev                = [];
        cfg_Q_ev.smooth         = 'gauss';
        cfg_Q_ev.gausswin_sd    = cfg.Q_gausswin_sd;
        cfg_Q_ev.dt             = ev_dt;
        cfg_Q_ev.tvec_edges     = tvec_ev_cont(1):ev_dt:tvec_ev_cont(end) + ev_dt;
        F_ev_full               = MakeQfromS(cfg_Q_ev, myCell);
        F_ev_full.data          = F_ev_full.data ./ cfg_Q_ev.dt;

        % Align FR to EV uniform grid
        F_idx_ev    = nearest_idx3(tvec_ev_cont, F_ev_full.tvec);
        F_idx_ev    = min(F_idx_ev, size(F_ev_full.data, 2));
        EV_FR_cont  = F_ev_full.data(1, F_idx_ev);

        % Non-NaN mask for scatter
        ev_valid_mask   = ~isnan(ev_cont_data);
        EV_scatter      = ev_cont_data(ev_valid_mask);
        EV_scatter_FR   = EV_FR_cont(ev_valid_mask);

        %% ============================================================
        %  STEP 5: Desaccade spike train for TuningCurves (separate path)
        %  ============================================================
        [~, keepSpk] = restrict(myCell, pre_sacc, post_sacc);
        myCellr      = myCell;
        myCellr.t{1,1}(keepSpk{1,1} == 1) = [];

        if cfg.remove_stationary
            myCellr = antirestrict(myCellr, stat_data.STtstart, stat_data.STtend);
        end

        %% ============================================================
        %  STEP 6: Compute tuning curves (spike train + tuning variable)
        %  ============================================================
        cfg_tcAHV           = [];
        cfg_tcAHV.nBins     = length(cfg.AHV_binEdges{1}) - 1;
        cfg_tcAHV.binEdges  = cfg.AHV_binEdges;
        cfg_tcAHV.minOcc    = cfg.AHV_minOcc;
        cfg_tcAHV.occ_dt    = occ_dt_AHV;
        tc_AHV              = TuningCurves(cfg_tcAHV, myCellr, AHVr);

        cfg_tcEV            = [];
        cfg_tcEV.nBins      = length(cfg.EV_binEdges{1}) - 1;
        cfg_tcEV.binEdges   = cfg.EV_binEdges;
        cfg_tcEV.minOcc     = cfg.EV_minOcc;
        cfg_tcEV.occ_dt     = ev_dt;
        tc_EV               = TuningCurves(cfg_tcEV, myCellr, diffHr_slow);

        %% ============================================================
        %  STEP 7: Panel A scatter data + regression
        %  ============================================================
        ev_idx_for_ahv  = nearest_idx3(AHVr.tvec, tvec_ev_cont);
        ev_idx_for_ahv  = min(ev_idx_for_ahv, length(ev_cont_data));
        AHV_for_scatter = AHVr.data(:)';
        EV_for_scatter  = ev_cont_data(ev_idx_for_ahv);

        valid_panelA    = ~isnan(AHV_for_scatter) & ~isnan(EV_for_scatter);
        ahv_valid       = AHV_for_scatter(valid_panelA)';
        ev_valid        = EV_for_scatter(valid_panelA)';
        [r_val, p_val]  = corr(ahv_valid, ev_valid);

        if cfg.robustFit
            lm          = fitlm(ahv_valid, ev_valid, 'RobustOpts', 'on');
            reg_slope   = lm.Coefficients.Estimate(2);
            reg_intcpt  = lm.Coefficients.Estimate(1);
            reg_method  = 'fitlm_robust';
            reg_R2      = lm.Rsquared.Ordinary;
        else
            X_design        = [ones(length(ahv_valid), 1), ahv_valid];
            [b_reg, ~, ~, ~, stats_reg] = regress(ev_valid, X_design);
            reg_slope       = b_reg(2);
            reg_intcpt      = b_reg(1);
            reg_method      = 'regress_OLS';
            reg_R2          = stats_reg(1);
        end

        if cfg.verbose
            fprintf('  Regression (%s): slope=%.4f  intercept=%.4f  R2=%.4f  r=%.4f  p=%.2e\n', ...
                reg_method, reg_slope, reg_intcpt, reg_R2, r_val, p_val);
        end

        %% Cache neuron data
        nd                  = [];
        nd.tfile_path       = tfile_path;
        nd.tfile_name       = tfile_name;
        nd.SSN              = SSN;
        nd.cell_idx         = cell_idx;
        nd.AHVr             = AHVr;
        nd.AHV_FR           = AHV_FR;
        nd.EV_scatter       = EV_scatter;
        nd.EV_scatter_FR    = EV_scatter_FR;
        nd.tc_AHV           = tc_AHV;
        nd.tc_EV            = tc_EV;
        nd.AHV_for_scatter  = AHV_for_scatter;
        nd.EV_for_scatter   = EV_for_scatter;
        nd.r_val            = r_val;
        nd.p_val            = p_val;
        nd.reg_slope        = reg_slope;
        nd.reg_intcpt       = reg_intcpt;
        nd.reg_method       = reg_method;
        nd.reg_R2           = reg_R2;
        nd.ev_params        = ev_params;
        nd.ev_out           = ev_out;
        nd.remove_stationary = cfg.remove_stationary;
        neuron_data{ii}     = nd;

        popdir;
    end % neuron loop
end % cache / compute

%% Plot all neurons
fig_handles = gobjects(nToPlot, 1);

for ii = 1:nToPlot
    nd = neuron_data{ii};
    if isempty(nd); continue; end

    %% Create figure
    if cfg.doShow
        fig_h = figure('Position', [50 200 cfg.figWidth cfg.figHeight]);
    else
        fig_h = figure('Position', [50 200 cfg.figWidth cfg.figHeight], ...
                       'Visible', 'off');
    end
    fig_handles(ii) = fig_h;

    % --- PANEL A: AHV vs EV scatterplot ---
    ax_A = subplot(1, 3, 1); hold on;

    if cfg.panelA_useGappy && isfield(nd, 'ev_out')
        % Use gappy (non-interpolated) EV — breaks up horizontal banding
        gappy_tvec = nd.ev_out.ev_gappy.tvec;
        gappy_data = nd.ev_out.ev_gappy.data;
        if size(gappy_tvec, 1) == 1; gappy_tvec = gappy_tvec'; end
        if size(gappy_data, 1) ~= 1; gappy_data = gappy_data'; end
        gappy_ahv_idx = nearest_idx3(gappy_tvec, nd.AHVr.tvec);
        gappy_ahv_idx = min(gappy_ahv_idx, length(nd.AHVr.data));
        panelA_ahv = nd.AHVr.data(gappy_ahv_idx);
        panelA_ev  = gappy_data;
    else
        panelA_ahv = nd.AHV_for_scatter;
        panelA_ev  = nd.EV_for_scatter;
    end

    plot(panelA_ahv, panelA_ev, '.', ...
        'MarkerSize', cfg.MarkerSize, 'Color', cfg.greyColor);

    % Regression line spans actual data range (not axis limits)
    valid_panelA_mask = ~isnan(nd.AHV_for_scatter) & ~isnan(nd.EV_for_scatter);
    ahv_data_range    = [min(nd.AHV_for_scatter(valid_panelA_mask)), ...
                         max(nd.AHV_for_scatter(valid_panelA_mask))];
    ahv_line_vals     = linspace(ahv_data_range(1), ahv_data_range(2), 100);
    ev_line_vals      = nd.reg_slope * ahv_line_vals + nd.reg_intcpt;
    plot(ahv_line_vals, ev_line_vals, 'k', 'LineWidth', 3);

    axis([cfg.panelA_AHV_lim cfg.panelA_EV_lim]);
    line([0 0], cfg.panelA_EV_lim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
    line(cfg.panelA_AHV_lim, [0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

    xlabel('AHV (deg/s)',          'FontSize', cfg.FontSize);
    ylabel('Eye velocity (deg/s)', 'FontSize', cfg.FontSize);
    set(ax_A, 'FontSize', cfg.FontSize - 2, 'TickDir', 'out', 'Box', 'off');

    text(0.05, 0.15, sprintf('r = %.2f', nd.r_val), ...
        'Units', 'normalized', 'FontSize', cfg.FontSize - 2, 'FontWeight', 'bold');
    if nd.p_val < 0.001
        p_str = 'p < 0.001';
    else
        p_str = sprintf('p = %.3f', nd.p_val);
    end
    text(0.05, 0.05, p_str, ...
        'Units', 'normalized', 'FontSize', cfg.FontSize - 2);

    % --- PANEL B: AHV tuning curve ---
    ax_B = subplot(1, 3, 2); hold on;
    plot(nd.AHVr.data, nd.AHV_FR, '.', ...
        'MarkerSize', cfg.MarkerSize, 'Color', cfg.greyColor);

    ahv_tc_vals = nd.tc_AHV.tc;
    ahv_good    = nd.tc_AHV.occ_hist > cfg.AHV_minOcc;
    if cfg.smoothTC
        ahv_tc_plot = smoothdata(ahv_tc_vals(ahv_good), 'gaussian', cfg.smoothTC_win);
    else
        ahv_tc_plot = ahv_tc_vals(ahv_good);
    end
    plot(nd.tc_AHV.binCenters(ahv_good), ahv_tc_plot, 'k', 'LineWidth', cfg.LineWidth);

    xlabel('AHV (deg/s)', 'FontSize', cfg.FontSize);
    set(ax_B, 'FontSize', cfg.FontSize - 2, 'TickDir', 'out', 'Box', 'off');
    set(ax_B, 'YTick', [], 'YTickLabel', []);  % no yticks — shared ylabel on Panel C
    text(0.10, 0.92, 'CW',  'Units', 'normalized', 'FontSize', cfg.FontSize, 'FontWeight', 'bold');
    text(0.78, 0.92, 'CCW', 'Units', 'normalized', 'FontSize', cfg.FontSize, 'FontWeight', 'bold');

    % --- PANEL C: EV tuning curve ---
    ax_C = subplot(1, 3, 3); hold on;
    plot(nd.EV_scatter, nd.EV_scatter_FR, '.', ...
        'MarkerSize', cfg.MarkerSize, 'Color', cfg.greyColor);

    ev_tc_vals  = nd.tc_EV.tc;
    ev_good     = nd.tc_EV.occ_hist > cfg.EV_minOcc;
    if cfg.smoothTC
        ev_tc_plot = smoothdata(ev_tc_vals(ev_good), 'gaussian', cfg.smoothTC_win);
    else
        ev_tc_plot = ev_tc_vals(ev_good);
    end
    plot(nd.tc_EV.binCenters(ev_good), ev_tc_plot, 'k', 'LineWidth', cfg.LineWidth);

    xlabel('Eye velocity (deg/s)', 'FontSize', cfg.FontSize);
    set(ax_C, 'FontSize', cfg.FontSize - 2, 'TickDir', 'out', 'Box', 'off');
    % Move shared FR ylabel to right side of Panel C
    set(ax_C, 'YAxisLocation', 'right');
    ylabel('FR (Hz)', 'FontSize', cfg.FontSize);
    text(0.10, 0.92, 'CW',  'Units', 'normalized', 'FontSize', cfg.FontSize, 'FontWeight', 'bold');
    text(0.78, 0.92, 'CCW', 'Units', 'normalized', 'FontSize', cfg.FontSize, 'FontWeight', 'bold');

    %% Match FR y-axes: Panels B and C share the same FR range
    all_FR_vals = [nd.AHV_FR(:); nd.EV_scatter_FR(:)];
    fr_max      = max(all_FR_vals(isfinite(all_FR_vals)));
    fr_max      = ceil(fr_max * 1.05);  % 5% headroom
    fr_ylim     = [0 fr_max];

    % Set axis limits (x and y) — do this once, then draw zero lines
    axis(ax_B, [cfg.panelA_AHV_lim fr_ylim]);
    axis(ax_C, [cfg.panelC_EV_lim  fr_ylim]);

    % Zero lines drawn after final axis limits are set
    axes(ax_B);
    line([0 0], fr_ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
    axes(ax_C);
    line([0 0], fr_ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

    %% Align panel tops and bottoms, and widen panels
    % Use Panel A's vertical position as reference, apply to B and C
    pos_A = get(ax_A, 'Position');
    pos_B = get(ax_B, 'Position');
    pos_C = get(ax_C, 'Position');
    % Widen each panel around its center
    for ax_pos = {'pos_A', 'pos_B', 'pos_C'}
        p = eval(ax_pos{1});
        old_w = p(3);
        new_w = old_w * cfg.panelWidthScale;
        p(1) = p(1) - (new_w - old_w) / 2;  % shift left to keep centered
        p(3) = new_w;
        eval([ax_pos{1} ' = p;']);
    end
    pos_B(2) = pos_A(2);          % match bottom
    pos_B(4) = pos_A(4);          % match height
    pos_C(2) = pos_A(2);
    pos_C(4) = pos_A(4);
    set(ax_A, 'Position', pos_A);
    set(ax_B, 'Position', pos_B);
    set(ax_C, 'Position', pos_C);

    %% Embed metadata
    fig_h.UserData.created_by       = 'plot_AHV_vs_EV_figure';
    fig_h.UserData.version          = '2.3';
    fig_h.UserData.date_created     = datestr(now, 'yyyy-mm-dd');
    fig_h.UserData.neuron_tfile     = nd.tfile_path;
    fig_h.UserData.session          = nd.SSN;
    fig_h.UserData.cell_idx         = nd.cell_idx;
    fig_h.UserData.r_value          = nd.r_val;
    fig_h.UserData.p_value          = nd.p_val;
    fig_h.UserData.regression_slope = nd.reg_slope;
    fig_h.UserData.regression_intcpt = nd.reg_intcpt;
    fig_h.UserData.regression_method = nd.reg_method;
    fig_h.UserData.regression_R2    = nd.reg_R2;
    fig_h.UserData.EV_pipeline      = 'compute_slow_phase_EV (Order B with stationary removal)';

    %% Save
    if cfg.robustFit
        fit_suffix = '_robust';
    else
        fit_suffix = '_OLS';
    end
    if cfg.remove_stationary
        stat_suffix = '_noStat';
    else
        stat_suffix = '';
    end
    fname_base = sprintf('figS2_%s%s%s', strrep(nd.tfile_name, '_', '-'), fit_suffix, stat_suffix);

    if cfg.doSavePNG
        print(fig_h, fullfile(save_subfolder, [fname_base '.png']), '-dpng', '-r300');
        if cfg.verbose
            fprintf('Saved PNG: %s.png\n', fname_base);
        end
    end

    if cfg.doSaveFIG
        savefig(fig_h, fullfile(save_subfolder, [fname_base '.fig']));
        if cfg.verbose
            fprintf('Saved FIG: %s.fig\n', fname_base);
        end
    end

    if cfg.doSaveEPS
        print(fig_h, fullfile(save_subfolder, [fname_base '.eps']), '-depsc', '-r300');
        if cfg.verbose
            fprintf('Saved EPS: %s.eps\n', fname_base);
        end
    end

    if ~cfg.doShow
        close(fig_h);
    end
end % plot loop

%% Save params and neuron data cache
if cfg.doSavePNG || cfg.doSaveFIG || cfg.doSaveEPS
    if cfg.robustFit
        mat_suffix = '_robust';
    else
        mat_suffix = '_OLS';
    end
    if cfg.remove_stationary
        mat_stat_suffix = '_noStat';
    else
        mat_stat_suffix = '';
    end
    mat_full_suffix = [mat_suffix mat_stat_suffix];
    save(fullfile(save_subfolder, ['figS2_params' mat_full_suffix '.mat']), 'params');
    save(fullfile(save_subfolder, ['figS2_neuron_data' mat_full_suffix '.mat']), ...
        'neuron_data', 'plot_indices', '-v7.3');
    if cfg.verbose
        fprintf('\nParams saved to: %s\n', save_subfolder);
        fprintf('Neuron data cached. To replot without recomputing:\n');
        fprintf('  cfg.cache_path = ''%s'';\n', ...
            fullfile(save_subfolder, ['figS2_neuron_data' mat_full_suffix '.mat']));
    end
end

end

%% Version History
% 2026-03-18  v2.3
%                   - FR computation refactored: MakeQfromS on RAW spike
%                     train first, then restrict FR vector with same masks
%                     as AHV/EV. Fixes 5000+ Hz artifact in noStat version.
%                   - TuningCurves still uses desaccaded spike train
%                     (separate path from FR scatter computation)
%                   - Cleaned up duplicate code blocks from v2.2
% 2026-03-18  v2.2
%                   - Added cfg.smoothTC_win (default 11, was hardcoded 5)
%                     for adjustable TC smoothing in workspace mode
%                   - Added cfg.remove_stationary (default 0). When 1,
%                     removes platform-stationary periods from both AHV
%                     (via antirestrict) and EV (passed to compute_slow_phase_EV).
%                     Requires recompute (not workspace mode).
%                   - Filenames include _noStat suffix when remove_stationary=1
% 2026-03-18  v2.1
%                   Cosmetic polish for publication:
%                   - panelA_EV_lim: [-200 200] -> [-50 50] (data is ±33)
%                   - panelA_AHV_lim: [-200 200] -> [-175 175]
%                   - AHV_minOcc: 2 -> 1 s (recover more edge bins in TC)
%                   - FontSize: 28 -> 32
%                   - MarkerSize: 4 -> 5
%                   - LineWidth: 5 -> 4 (slightly thinner TC for cleaner look)
%                   - greyColor: [0.45 0.45 0.45] -> [0.5 0.5 0.5] (lighter)
%                   - Added panel labels (A, B, C) in upper-left
%                   - Regression line now spans actual data range, not ±200
%                   - figHeight: 550 -> 600 (more vertical room)
%                   - Zero-line drawing consolidated (drawn once after ylim)
%                   - AHV bin edges now match panelA_AHV_lim (was ±200)
%                   - Added workspace mode: cfg.neuron_data skips recompute
%                   - 3rd output returns neuron_data for workspace reuse
%                   - Regression: polyfit -> regress (OLS) with cfg.robustFit
%                     option for fitlm robust (bisquare). Stores slope,
%                     intercept, R2, method in neuron_data.
% 2026-03-17  v2.0
%                   - EV pipeline now via compute_slow_phase_EV (Order B
%                     with stationary-period removal). Fixes the missing
%                     stationary exclusion that left many near-zero EV
%                     points in the scatter and TC.
%                   - AHV also now has stationary periods removed
%                     (antirestrict with STtstart/STtend).
%                   - minOcc defaults changed from sample counts to seconds
%                     (AHV_minOcc: 100 → 2, EV_minOcc: 25 → 0.5) to match
%                     TuningCurves occ_hist units.
%                   - Removed cfg fields now handled by compute_slow_phase_EV:
%                     EV_velocity_threshold, medfilt_win, gap_thresh
% 2026-03-17  v1.3
%                   - Fixed nearest_idx3 argument order in EV FR alignment
% 2026-03-17  v1.2
%                   - Panel C: interpolate onto uniform grid, NaN across gaps
%                   - Panel A scatter uses interpolated EV
%                   - Cache stores EV_scatter/EV_scatter_FR
% 2026-03-17  v1.1
%                   - EV pipeline changed to Order B
%                   - Panels B and C share matched FR y-axis
%                   - Added varargin support
% 2026-03-17  v1.0 Initial version.
