function [fig_handle, params, cell_data] = plot_fig4B_ver_1_5(varargin)
% plot_fig4B_ver_1_5  New Figure 4B: EV tuning curve + preferred saccade PETH.
%
%   Generates a 3-row x 2-column figure panel for the manuscript:
%     Row 1: "Same" example neuron (slow-phase and quick-phase same direction)
%     Row 2: "Opposite" example neuron (slow-phase and quick-phase opposite)
%     Row 3: "Bidirectional" example neuron (slow-phase linear, quick-phase v-shape)
%
%   Left column:  Slow-phase EV tuning curve (grey scatter + black TC line)
%                  CW / CCW direction labels, preferred direction bolded.
%   Right column: Preferred saccade raster (top) + smoothed FR PETH (bottom)
%                  Title indicates saccade direction shown.
%
%   The visual comparison between the bolded EV direction and the saccade
%   PETH direction title makes same/opposite/bidirectional immediately clear.
%
%   INPUTS:
%       varargin    - name-value pairs (see defaults below)
%
%   KEY OPTIONS:
%       neuron_same         - neuron name for Same row
%       neuron_opposite     - neuron name for Opposite row
%       neuron_bidirectional - neuron name for Bidirectional row
%       cell_data           - pre-computed cell data (skip reloading)
%       tfilelist_path      - path to NPH_all_tfilelist.mat
%       data_root           - root data path
%       save_root           - analysis output path
%       doSavePNG, doSaveFIG, doShow
%       FontSize            - base font size (default: 30)
%
%   OUTPUTS:
%       fig_handle  - figure handle
%       params      - struct of all parameter values used
%       cell_data   - cell array {3x1} of per-neuron data (pass back to skip recompute)
%
%   USAGE:
%       % First run (from scratch):
%       [fh, par, cd] = plot_fig4B_ver_1_5();
%
%       % Subsequent runs (workspace mode):
%       [fh2, par2] = plot_fig4B_ver_1_5('cell_data', cd, 'doShow', 1);
%
%   DEPENDENCIES:
%       compute_slow_phase_EV, LoadSpikesQuick, SelectTS, TuningCurves,
%       MakeQfromS, nearest_idx3, tsd, restrict, pushdir, popdir,
%       HD_GetSSN, process_varargin, ProcessConfig2
%

%% ========================================================================
%  DEFAULTS
%  ========================================================================
neuron_same             = 'M079-2020-10-28-TT07_2';
neuron_opposite         = 'M085-2020-11-13-TT01_1';
neuron_bidirectional    = 'M086-2020-11-20-TT02_1';

tfilelist_path  = 'C:\Jeff\U01\datatouse\NPH all tfilelist.mat';
data_root       = 'C:\Jeff\U01\datatouse';
save_root       = 'C:\Jeff\U01\datatouse\___analysis_output';

% EV TC parameters
EV_binEdges         = linspace(-40, 40, 81);
EV_minOcc           = 0.5;          % seconds
smoothTC_win        = 11;           % bins, Gaussian smoothdata window
velocity_threshold  = 35;           % deg/s
Q_gausswin_sd       = 0.05;         % s, Gaussian FR smoothing kernel SD
saccade_pre         = 0.1;          % s to exclude before saccade
saccade_post        = 0.1;          % s to exclude after saccade
remove_stationary   = 0;            % 0 = keep stationary periods

% Saccade PETH parameters
peth_window         = [-0.2 0.1];   % s, matches manuscript methods
peth_dt             = 0.005;        % s, bin width
peth_smooth_win     = 5;            % bins, movmean window before smoothdata
peth_excessBounds   = 1;            % s, extra restrict margin (matches saccadePETH_one_neuron)

% Figure parameters
FontSize            = 30;
MarkerSize_scatter  = 4;            % grey EV scatter dots
MarkerSize_raster   = 3;            % raster tick marks
LineWidth_TC        = 5;            % EV TC overlay line
LineWidth_FR        = 4;            % FR PETH line
greyColor           = [0.5 0.5 0.5];
rasterColor         = [0.15 0.15 0.15];
fr_color_same       = [0.30 0.69 0.29];     % green (matches Fig4A green rows)
fr_color_opposite   = [0.90 0.45 0.15];     % orange (matches Fig4A orange rows)
fr_color_bidir      = [0.55 0.70 0.85];     % light blue (matches Fig4A blue rows)
figWidth            = 1600;
figHeight           = 850;

% Output control
doSavePNG   = 1;
doSaveFIG   = 0;
doShow      = 0;
verbose     = 1;

% Cache
cell_data   = {};

%% Process varargin
process_varargin(varargin);

%% ========================================================================
%  PARAMS
%  ========================================================================
params = struct();
params.function_name        = 'plot_fig4B';
params.version              = '1.5';
params.date_run             = datestr(now, 'yyyy-mm-dd');
params.neuron_same          = neuron_same;
params.neuron_opposite      = neuron_opposite;
params.neuron_bidirectional = neuron_bidirectional;
params.EV_binEdges          = EV_binEdges;
params.EV_minOcc            = EV_minOcc;
params.smoothTC_win         = smoothTC_win;
params.velocity_threshold   = velocity_threshold;
params.Q_gausswin_sd        = Q_gausswin_sd;
params.saccade_pre          = saccade_pre;
params.saccade_post         = saccade_post;
params.remove_stationary    = remove_stationary;
params.peth_window          = peth_window;
params.peth_dt              = peth_dt;
params.peth_smooth_win      = peth_smooth_win;
params.peth_excessBounds    = peth_excessBounds;
params.FontSize             = FontSize;

%% ========================================================================
%  NEURON INFO TABLE (from Chenhao's fig4_table)
%  Hardcoded category assignments for the 3 default neurons.
%  If user provides different neurons, they must also provide category info
%  or this function reads it from the table.
%  ========================================================================
% Category info: {name, slow_pref ('CCW'|'CW'), quick_pref ('CCW'|'CW'|'Burst'), row_label}
neuron_list = {
    neuron_same,            'Same';
    neuron_opposite,        'Opposite';
    neuron_bidirectional,   'Bidirectional';
};
nNeurons = size(neuron_list, 1);

% Direction assignments for default neurons (from Chenhao's table)
% Map: neuron_name -> {slow_pref, quick_pref}
dir_map = containers.Map();
dir_map('M079-2020-10-28-TT07_2')  = struct('slow_pref', 'CCW', 'quick_pref', 'CCW');
dir_map('M085-2020-11-13-TT01_1')  = struct('slow_pref', 'CCW', 'quick_pref', 'CW');
dir_map('M086-2020-11-20-TT02_1')  = struct('slow_pref', 'CCW', 'quick_pref', 'Burst');
dir_map('M052-2020-09-04-TT03_1')  = struct('slow_pref', 'CW',  'quick_pref', 'CW');
dir_map('M222-2022-01-21-TT01_1')  = struct('slow_pref', 'CW',  'quick_pref', 'CW');
dir_map('M039-2020-08-21-1-TT01_2') = struct('slow_pref', 'CCW', 'quick_pref', 'CW');
dir_map('M039-2020-08-21-1-TT01_1') = struct('slow_pref', 'CCW', 'quick_pref', 'Burst');
dir_map('M090-2020-10-02-TT05_1')  = struct('slow_pref', 'CCW', 'quick_pref', 'CCW');
dir_map('M287-2022-08-09-1-TT07_1') = struct('slow_pref', 'CCW', 'quick_pref', 'CCW');
dir_map('M079-2020-10-28-TT07_1')  = struct('slow_pref', 'CCW', 'quick_pref', 'CW');
dir_map('M085-2020-11-13-TT08_1')  = struct('slow_pref', 'CW',  'quick_pref', 'CCW');
dir_map('M247-2022-01-10-TT01_1')  = struct('slow_pref', 'CW',  'quick_pref', 'CCW');
dir_map('M090-2020-10-05-TT03_1')  = struct('slow_pref', 'CW',  'quick_pref', 'Burst');
dir_map('M247-2022-01-14-1-TT07_1') = struct('slow_pref', 'CW', 'quick_pref', 'Burst');

%% ========================================================================
%  LOAD TFILELIST & COMPUTE DATA (or use cache)
%  ========================================================================
if ~isempty(cell_data) && length(cell_data) == nNeurons
    if verbose; fprintf('Using pre-computed cell_data from workspace.\n'); end
else
    %% Load master tfilelist
    tfl = load(tfilelist_path, 'NPH_all_tfilelist');
    NPH_all_tfilelist = tfl.NPH_all_tfilelist;

    cell_data = cell(nNeurons, 1);

    for iN = 1:nNeurons
        neuron_name = neuron_list{iN, 1};
        row_label   = neuron_list{iN, 2};

        if verbose
            fprintf('\n=== Neuron %d/%d: %s (%s) ===\n', iN, nNeurons, neuron_name, row_label);
        end

        %% Find session path from tfilelist
        sess_path = '';
        tfile_full_path = '';
        for kk = 1:length(NPH_all_tfilelist)
            [kk_folder, kk_name] = fileparts(NPH_all_tfilelist{kk});
            if strcmp(kk_name, neuron_name)
                sess_path = kk_folder;
                tfile_full_path = NPH_all_tfilelist{kk};
                break
            end
        end
        if isempty(sess_path)
            error('Could not find neuron %s in tfilelist.', neuron_name);
        end

        %% cd to session and load spikes
        pushdir(sess_path);
        SSN = HD_GetSSN;
        [S_all, spike_cell_names, spk_params] = LoadSpikesQuick();

        % Find cell index (LoadSpikesQuick preserves underscores)
        cell_idx = [];
        for cc = 1:length(spike_cell_names)
            if strcmp(spike_cell_names{cc}, neuron_name)
                cell_idx = cc;
                break
            end
        end
        if isempty(cell_idx)
            % Try matching without session prefix
            % e.g., spike_cell_names might be full SSN-TT format
            for cc = 1:length(spike_cell_names)
                if contains(spike_cell_names{cc}, strrep(neuron_name, [SSN '-'], ''))
                    cell_idx = cc;
                    break
                end
            end
        end
        if isempty(cell_idx)
            error('Could not match %s to loaded spike trains. Names: %s', ...
                neuron_name, strjoin(spike_cell_names, ', '));
        end
        myCell_raw = SelectTS([], S_all, cell_idx);

        %% Load saccade data
        sacc_file = strcat(SSN, '-saccades-edited.mat');
        sacc_data = load(sacc_file, 'combinedSaccades', 'temporalSaccades', ...
            'nasalSaccades', 'diffHdeg');
        temporal_sacc_times = sacc_data.temporalSaccades;
        nasal_sacc_times    = sacc_data.nasalSaccades;
        combinedSacc        = sacc_data.combinedSaccades;

        %% ==============================================================
        %  EV TUNING CURVE COMPUTATION
        %  ==============================================================

        % Step A: Compute slow-phase EV
        [ev_out, ~] = compute_slow_phase_EV('sess_dir', sess_path, ...
            'saccade_pre', saccade_pre, 'saccade_post', saccade_post, ...
            'velocity_threshold', velocity_threshold, ...
            'remove_stationary', remove_stationary, ...
            'verbose', 0);

        tvec_ev_interp  = ev_out.ev_interp.tvec;
        data_ev_interp  = ev_out.ev_interp.data;
        ev_dt_val       = ev_out.ev_dt;

        % Ensure correct orientation
        if size(tvec_ev_interp, 1) == 1; tvec_ev_interp = tvec_ev_interp'; end
        if size(data_ev_interp, 1) ~= 1; data_ev_interp = data_ev_interp'; end

        % Step B: Compute FR from RAW spike train on EV timebase
        cfg_Q_ev                = [];
        cfg_Q_ev.smooth         = 'gauss';
        cfg_Q_ev.gausswin_sd    = Q_gausswin_sd;
        cfg_Q_ev.dt             = ev_dt_val;
        cfg_Q_ev.tvec_edges     = tvec_ev_interp(1) : ev_dt_val : tvec_ev_interp(end) + ev_dt_val;
        F_ev                    = MakeQfromS(cfg_Q_ev, myCell_raw);
        F_ev.data               = F_ev.data ./ cfg_Q_ev.dt;

        % Align FR to EV uniform grid
        F_idx_ev    = nearest_idx3(tvec_ev_interp, F_ev.tvec);
        F_idx_ev    = min(F_idx_ev, size(F_ev.data, 2));
        FR_ev_aligned = F_ev.data(1, F_idx_ev);

        % Non-NaN mask for scatter
        ev_valid        = ~isnan(data_ev_interp);
        ev_scatter_vals = data_ev_interp(ev_valid);
        fr_scatter_vals = FR_ev_aligned(ev_valid);

        % Step C: Desaccade spike train for TuningCurves (separate path)
        sacc_times_clean = combinedSacc(~isnan(combinedSacc));
        pre_sacc_win     = sacc_times_clean - saccade_pre;
        post_sacc_win    = sacc_times_clean + saccade_post;
        [~, keepSpk]     = restrict(myCell_raw, pre_sacc_win, post_sacc_win);
        myCell_desacc    = myCell_raw;
        myCell_desacc.t{1,1}(keepSpk{1,1} == 1) = [];

        % Step D: Compute EV tuning curve
        diffHr_slow_tsd = tsd(tvec_ev_interp, data_ev_interp);

        cfg_tcEV            = [];
        cfg_tcEV.nBins      = length(EV_binEdges) - 1;
        cfg_tcEV.binEdges   = {EV_binEdges};
        cfg_tcEV.minOcc     = EV_minOcc;
        cfg_tcEV.occ_dt     = ev_dt_val;
        tc_EV_result        = TuningCurves(cfg_tcEV, myCell_desacc, diffHr_slow_tsd);

        % Smooth TC for well-occupied bins
        good_bins_ev        = tc_EV_result.occ_hist > EV_minOcc;
        tc_ev_smooth_vals   = NaN(1, length(tc_EV_result.tc));
        tc_ev_smooth_vals(good_bins_ev) = smoothdata(tc_EV_result.tc(good_bins_ev), ...
            'gaussian', smoothTC_win);

        if verbose
            fprintf('  EV TC: %d good bins out of %d\n', sum(good_bins_ev), length(good_bins_ev));
        end

        %% ==============================================================
        %  SACCADE PETH COMPUTATION (inline, both directions)
        %  ==============================================================
        peth_bins    = linspace(peth_window(1), peth_window(2), diff(peth_window)/peth_dt + 1);
        peth_centers = peth_bins(1:end-1) + diff(peth_bins) / 2;

        % Temporal saccade PETH
        [spk_rel_T, evt_idx_T] = compute_raster_data(myCell_raw, ...
            temporal_sacc_times, peth_window, peth_excessBounds);
        [m_T, ~] = histcounts(spk_rel_T, peth_bins);
        nTemporalEvents = length(temporal_sacc_times);
        fr_peth_T = m_T / peth_dt / max(nTemporalEvents, 1);
        fr_peth_T_smooth = smoothdata(movmean(fr_peth_T, peth_smooth_win), 'gaussian', peth_smooth_win);

        % Nasal saccade PETH
        [spk_rel_N, evt_idx_N] = compute_raster_data(myCell_raw, ...
            nasal_sacc_times, peth_window, peth_excessBounds);
        [m_N, ~] = histcounts(spk_rel_N, peth_bins);
        nNasalEvents = length(nasal_sacc_times);
        fr_peth_N = m_N / peth_dt / max(nNasalEvents, 1);
        fr_peth_N_smooth = smoothdata(movmean(fr_peth_N, peth_smooth_win), 'gaussian', peth_smooth_win);

        %% Determine preferred directions
        if dir_map.isKey(neuron_name)
            dir_info = dir_map(neuron_name);
        else
            error('No direction info for neuron %s. Add to dir_map or provide fig4_table.', ...
                neuron_name);
        end

        % Map quick phase preference to temporal/nasal
        % CCW saccade = temporal saccade, CW saccade = nasal saccade
        if strcmp(dir_info.quick_pref, 'CCW')
            pref_sacc_label   = 'CCW';
            pref_spk_rel      = spk_rel_T;
            pref_evt_idx      = evt_idx_T;
            pref_fr_smooth    = fr_peth_T_smooth;
            pref_nEvents      = nTemporalEvents;
        elseif strcmp(dir_info.quick_pref, 'CW')
            pref_sacc_label   = 'CW';
            pref_spk_rel      = spk_rel_N;
            pref_evt_idx      = evt_idx_N;
            pref_fr_smooth    = fr_peth_N_smooth;
            pref_nEvents      = nNasalEvents;
        else  % Burst — pick direction matching slow-phase preference
            if strcmp(dir_info.slow_pref, 'CCW')
                pref_sacc_label   = 'CCW';
                pref_spk_rel      = spk_rel_T;
                pref_evt_idx      = evt_idx_T;
                pref_fr_smooth    = fr_peth_T_smooth;
                pref_nEvents      = nTemporalEvents;
            else
                pref_sacc_label   = 'CW';
                pref_spk_rel      = spk_rel_N;
                pref_evt_idx      = evt_idx_N;
                pref_fr_smooth    = fr_peth_N_smooth;
                pref_nEvents      = nNasalEvents;
            end
        end

        %% Cache all computed data
        nd = struct();
        nd.neuron_name      = neuron_name;
        nd.row_label        = row_label;
        nd.SSN              = SSN;
        nd.cell_idx         = cell_idx;
        nd.tfile_path       = tfile_full_path;
        nd.dir_info         = dir_info;
        % EV TC data
        nd.ev_scatter_vals  = ev_scatter_vals;
        nd.fr_scatter_vals  = fr_scatter_vals;
        nd.tc_binCenters    = tc_EV_result.binCenters;
        nd.tc_good_bins     = good_bins_ev;
        nd.tc_smooth_vals   = tc_ev_smooth_vals;
        nd.tc_occ_hist      = tc_EV_result.occ_hist;
        % Saccade PETH data
        nd.peth_centers         = peth_centers;
        nd.pref_sacc_label      = pref_sacc_label;
        nd.pref_spk_rel         = pref_spk_rel;
        nd.pref_evt_idx         = pref_evt_idx;
        nd.pref_fr_smooth       = pref_fr_smooth;
        nd.pref_nEvents         = pref_nEvents;
        % Also store both directions (for later flexibility)
        nd.fr_peth_T_smooth     = fr_peth_T_smooth;
        nd.fr_peth_N_smooth     = fr_peth_N_smooth;
        nd.spk_rel_T            = spk_rel_T;
        nd.evt_idx_T            = evt_idx_T;
        nd.spk_rel_N            = spk_rel_N;
        nd.evt_idx_N            = evt_idx_N;
        nd.nTemporalEvents      = nTemporalEvents;
        nd.nNasalEvents         = nNasalEvents;

        cell_data{iN} = nd;

        popdir;

        if verbose
            fprintf('  Preferred EV: %s, Preferred saccade: %s (%s)\n', ...
                dir_info.slow_pref, pref_sacc_label, row_label);
        end
    end
end

%% ========================================================================
%  CREATE FIGURE
%  ========================================================================
if doShow
    fig_handle = figure('Position', [50 50 figWidth figHeight]);
else
    fig_handle = figure('Position', [50 50 figWidth figHeight], 'Visible', 'off');
end

% Color per row
row_colors = {fr_color_same, fr_color_opposite, fr_color_bidir};

%% Manual axes positioning
left_margin     = 0.11;
right_margin    = 0.03;
top_margin      = 0.06;
bottom_margin   = 0.10;
row_gap         = 0.09;        % more breathing room between rows
col_gap         = 0.07;

total_height    = 1 - top_margin - bottom_margin;
row_height      = (total_height - (nNeurons - 1) * row_gap) / nNeurons;
total_width     = 1 - left_margin - right_margin;
col_w_left      = total_width * 0.42;
col_w_right     = total_width * 0.50;

ax_ev       = gobjects(nNeurons, 1);
ax_sacc     = gobjects(nNeurons, 1);    % single combined raster+FR axis

for iRow = 1:nNeurons
    row_bottom = bottom_margin + (nNeurons - iRow) * (row_height + row_gap);
    right_x    = left_margin + col_w_left + col_gap;

    % EV TC axes
    ax_ev(iRow) = axes('Position', [left_margin, row_bottom, col_w_left, row_height]);

    % Combined raster + FR axes (single panel)
    ax_sacc(iRow) = axes('Position', [right_x, row_bottom, col_w_right, row_height]);
end

%% ========================================================================
%  PLOT EACH ROW
%  ========================================================================
for iRow = 1:nNeurons
    nd = cell_data{iRow};

    %% ---- LEFT: EV Tuning Curve ----------------------------------------
    axes(ax_ev(iRow)); %#ok<LAXES>
    hold on;

    % Grey scatter
    plot(nd.ev_scatter_vals, nd.fr_scatter_vals, '.', ...
        'MarkerSize', MarkerSize_scatter, 'Color', greyColor);

    % Black TC line (only good bins)
    plot(nd.tc_binCenters(nd.tc_good_bins), nd.tc_smooth_vals(nd.tc_good_bins), ...
        'k', 'LineWidth', LineWidth_TC);

    % Vertical zero line
    yl_ev = ylim;
    line([0 0], yl_ev, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

    % Axis formatting — strip all ticks
    xlim([EV_binEdges(1) EV_binEdges(end)]);
    set(gca, 'FontSize', FontSize - 2, 'Box', 'off', ...
        'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');

    % CW / CCW labels with preferred bolded
    if strcmp(nd.dir_info.slow_pref, 'CW')
        cw_weight  = 'bold';
        ccw_weight = 'normal';
    else
        cw_weight  = 'normal';
        ccw_weight = 'bold';
    end
    text(0.06, 0.92, 'CW',  'Units', 'normalized', 'FontSize', FontSize, ...
        'FontWeight', cw_weight, 'Color', 'k');
    text(0.80, 0.92, 'CCW', 'Units', 'normalized', 'FontSize', FontSize, ...
        'FontWeight', ccw_weight, 'Color', 'k');

    % Row label (category) on left margin
    text(-0.22, 0.5, nd.row_label, 'Units', 'normalized', ...
        'FontSize', FontSize + 4, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'Rotation', 90, ...
        'Color', row_colors{iRow});

    % X-axis label only on bottom row
    if iRow == nNeurons
        set(gca, 'XTick', [-30 -15 0 15 30], 'FontSize', FontSize - 10);
        xlabel('Eye velocity (deg/s)', 'FontSize', FontSize - 4);
    end

    %% ---- RIGHT: Combined Raster + FR (yyaxis) -------------------------
    axes(ax_sacc(iRow)); %#ok<LAXES>

    % Left y-axis: raster dots
    yyaxis left;
    hold on;
    if ~isempty(nd.pref_spk_rel)
        plot(nd.pref_spk_rel, nd.pref_evt_idx, '.', ...
            'MarkerSize', MarkerSize_raster, 'Color', [0.55 0.55 0.55]);
    end
    ylim([0 nd.pref_nEvents + 1]);
    set(gca, 'YTick', [], 'YColor', 'none');

    % Vertical line at saccade onset (draw on left axis so it spans full height)
    line([0 0], [0 nd.pref_nEvents + 1], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);

    % Right y-axis: FR PETH line
    yyaxis right;
    hold on;
    plot(nd.peth_centers, nd.pref_fr_smooth, '-', ...
        'LineWidth', LineWidth_FR, 'Color', row_colors{iRow});
    fr_peak = max(nd.pref_fr_smooth);
    if fr_peak > 0
        ylim([0 fr_peak * 1.25]);
    end
    set(gca, 'YTick', [], 'YColor', 'none');

    % Common formatting
    xlim(peth_window);
    set(gca, 'Box', 'off', 'XTick', []);

    % X-axis label only on bottom row
    if iRow == nNeurons
        set(gca, 'XTick', [-0.2 -0.1 0 0.1], 'FontSize', FontSize - 10);
        xlabel('Time from saccade (s)', 'FontSize', FontSize - 4);
    end

    % Title: saccade direction (above the panel with clearance)
    title(sprintf('%s Saccade', nd.pref_sacc_label), ...
        'FontSize', FontSize, 'FontWeight', 'bold', 'Color', row_colors{iRow});
end

%% ========================================================================
%  METADATA
%  ========================================================================
fig_handle.UserData.created_by      = 'plot_fig4B_ver_1_5';
fig_handle.UserData.date_created    = datestr(now, 'yyyy-mm-dd');
fig_handle.UserData.params          = params;
fig_handle.UserData.neurons         = {cell_data{1}.neuron_name, ...
    cell_data{2}.neuron_name, cell_data{3}.neuron_name};

%% ========================================================================
%  SAVE
%  ========================================================================
if doSavePNG || doSaveFIG
    % Find next run number
    run_date   = datestr(now, 'yyyy-mm-dd');
    run_prefix = ['fig4B_' run_date '_run'];
    existing   = dir(fullfile(save_root, [run_prefix '*']));
    run_num    = length(existing) + 1;
    save_folder = fullfile(save_root, sprintf('%s%d', run_prefix, run_num));
    mkdir(save_folder);

    fname_base = 'fig4B_EV_TC_saccade_PETH';

    if doSavePNG
        print(fig_handle, fullfile(save_folder, [fname_base '.png']), '-dpng', '-r300');
        if verbose; fprintf('\nSaved PNG: %s\n', fullfile(save_folder, [fname_base '.png'])); end
    end
    if doSaveFIG
        savefig(fig_handle, fullfile(save_folder, [fname_base '.fig']));
        if verbose; fprintf('Saved FIG: %s\n', fullfile(save_folder, [fname_base '.fig'])); end
    end

    % Save params and cell_data
    save(fullfile(save_folder, 'fig4B_params.mat'), 'params', '-v7');
    save(fullfile(save_folder, 'fig4B_cell_data.mat'), 'cell_data', '-v7.3');

    if verbose
        fprintf('Data saved to: %s\n', save_folder);
        fprintf('To replot without recomputing:\n');
        fprintf('  [fh, par] = plot_fig4B_ver_1_5(''cell_data'', cell_data, ''doShow'', 1);\n');
    end
end

%% Print neuron IDs to command window
fprintf('\n--- Figure 4B Neuron Summary ---\n');
for iP = 1:nNeurons
    nd_print = cell_data{iP};
    fprintf('  Row %d (%s): %s  |  EV pref: %s  |  Saccade pref: %s\n', ...
        iP, nd_print.row_label, nd_print.neuron_name, ...
        nd_print.dir_info.slow_pref, nd_print.pref_sacc_label);
end
fprintf('--------------------------------\n');

if ~doShow
    close(fig_handle);
end

end

%% ========================================================================
%  LOCAL HELPER: compute raster data for one saccade direction
%  ========================================================================
function [spk_times_rel, event_indices] = compute_raster_data(myCell, saccade_times, window, excessBounds)
% Mimics saccadePETH_one_neuron's raster computation without plotting.
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
% 2026-04-02  v1.5  Fullscreen landscape figure (1600x850). Matched tick
%                    font sizes across both columns (FontSize-10). Xlabel
%                    font reduced to FontSize-4 for proportionality.
% 2026-04-02  v1.4  Resized figure to publication aspect ratio (750x900).
%                    Added xlabels on bottom row for both columns.
%                    Bottom margin increased for label clearance.
% 2026-04-02  v1.3  Combined raster and FR into single panel via yyaxis.
%                    Raster dots on left axis, FR line on right axis overlaid.
%                    More vertical space between rows. Removed separate FR subplot.
% 2026-04-02  v1.2  Stripped all ticks/tick labels for clean visual evaluation.
%                    Raster changed from dots to vectorized vertical tick marks.
%                    Saccade direction title repositioned above raster.
%                    Category labels enlarged. Removed xlabel/ylabel clutter.
% 2026-04-02  v1.1  Cleaned up labels: removed redundant ylabels from right
%                    panels and all but middle-row left panel. Removed neuron
%                    ID titles from figure, prints them to command window instead.
% 2026-04-02  v1.0  Initial version. Replaces AHV TC with EV TC in Figure 4B.
%                    Shows preferred saccade raster + FR PETH. Three example
%                    rows: same, opposite, bidirectional. EV TC via
%                    compute_slow_phase_EV (Order B). Saccade PETH computed
%                    inline. Category assignments from Chenhao's fig4_table.
