function [all_cell_data, pop_data, params] = plot_optotag_population_PETHs(varargin)
% plot_optotag_population_PETHs  v3.3
% Laser + saccade PETHs for 7 optotagged DTN-projecting NPH neurons.
%
%   2x3 layout:
%     Row 1: VOR Preferred | VOR Non-Preferred | Laser FR PETH
%     Row 2: Spont Preferred | Spont Non-Preferred | Laser Eye Position PETH
%
%   Colors match analyze_stationary_saccade_bursts_ver_1_3.
%   Laser FR panel: blue background shading, black FR line.
%   Laser eye position panel: baseline-subtracted tsdH aligned to laser onset.

%% =====================================================================
%  DEFAULTS
%  =====================================================================
precomputed       = [];
save_dir          = 'C:\Jeff\U01\datatouse\___analysis_output';
data_root         = 'C:\Jeff\U01\datatouse';

opto_neurons = {'M271-2021-08-30-TT02_1', 'M282-2022-02-10-TT03_2', ...
    'M282-2022-02-15-TT07_1', 'M391-2022-12-07-1-TT01_1', ...
    'M391-2022-12-07-2-TT02_1', 'M402-2022-11-02-2-TT02_3', ...
    'M444-2023-08-18-TT07_1'};
opto_sessions = {'M271-2021-08-30', 'M282-2022-02-10', 'M282-2022-02-15', ...
    'M391-2022-12-07-1', 'M391-2022-12-07-2', 'M402-2022-11-02-2', ...
    'M444-2023-08-18'};

laser_window      = [-1 2];
saccade_window    = [-0.2 0.2];
peth_dt           = 0.005;
peth_smooth_win   = 5;
peth_excessBounds = 1;
laser_baseline    = [-1 -0.2];
sacc_baseline     = [-0.2 -0.1];
burst_window      = [-0.05 0.05];
eye_baseline      = [-1 -0.2];   % for baseline-subtracting eye position

% Colors (matching v1.3 population figure)
color_vor_pref      = [0.20 0.50 0.80];
color_vor_nonpref   = [0.90 0.40 0.15];
color_spont_pref    = [0.08 0.25 0.50];
color_spont_nonpref = [0.60 0.22 0.05];
shade_vor_pref      = [0.72 0.84 0.96];
shade_vor_nonpref   = [0.97 0.80 0.65];
shade_spont_pref    = [0.30 0.45 0.65];
shade_spont_nonpref = [0.78 0.52 0.38];
color_laser_line    = [0 0 0];
color_laser_bg      = [0.4 0.7 1.0];
color_eye_line      = [0.60 0.15 0.60];   % purple for eye position
color_eye_shade     = [0.80 0.60 0.80];
shade_alpha         = 0.35;

FontSize = 20; TitleFontSize = 16; TickFontSize = 13; LineWidth_FR = 3;
doSavePNG = 1; doSaveFIG = 1; doSaveEPS = 1; doShow = 1;

process_varargin(varargin);
for iV = 1:2:length(varargin)
    if ischar(varargin{iV}) && strcmp(varargin{iV}, 'precomputed')
        precomputed = varargin{iV+1};
    end
end

%% =====================================================================
%  PARAMS
%  =====================================================================
params = struct();
params.function_name  = 'plot_optotag_population_PETHs';
params.version        = '3.3';
params.date_run       = datestr(now, 'yyyy-mm-dd HH:MM:SS');
params.laser_window   = laser_window;
params.saccade_window = saccade_window;
params.eye_baseline   = eye_baseline;
params.peth_dt        = peth_dt;
params.saccade_source = 'Chenhao onset times';

%% Output folder
run_date_str = datestr(now, 'yyyy-mm-dd');
run_num = 1;
out_folder = fullfile(save_dir, sprintf('optotag_pop_PETHs_%s_run%d', run_date_str, run_num));
while exist(out_folder, 'dir')
    run_num = run_num + 1;
    out_folder = fullfile(save_dir, sprintf('optotag_pop_PETHs_%s_run%d', run_date_str, run_num));
end
if doSavePNG || doSaveFIG || doSaveEPS; mkdir(out_folder); end
params.output_folder = out_folder;

%% Bin setup
laser_bins    = linspace(laser_window(1), laser_window(2), diff(laser_window)/peth_dt + 1);
laser_centers = laser_bins(1:end-1) + diff(laser_bins)/2;
nBins_laser   = length(laser_centers);

sacc_bins     = linspace(saccade_window(1), saccade_window(2), diff(saccade_window)/peth_dt + 1);
sacc_centers  = sacc_bins(1:end-1) + diff(sacc_bins)/2;
nBins_sacc    = length(sacc_centers);

% Eye position uses same time grid as laser PETH
eye_dt   = median(diff(laser_centers));
eye_grid = laser_centers;   % same grid
nBins_eye = nBins_laser;
eye_bl_mask = eye_grid >= eye_baseline(1) & eye_grid <= eye_baseline(2);

laser_bl_mask   = laser_centers >= laser_baseline(1) & laser_centers <= laser_baseline(2);
sacc_bl_mask    = sacc_centers >= sacc_baseline(1) & sacc_centers <= sacc_baseline(2);
sacc_burst_mask = sacc_centers >= burst_window(1) & sacc_centers <= burst_window(2);

nCells = length(opto_neurons);

%% =====================================================================
%  WORKSPACE MODE
%  =====================================================================
if ~isempty(precomputed)
    fprintf('=== WORKSPACE MODE (%d cells) ===\n', length(precomputed));
    all_cell_data = precomputed;
    pop_data = compute_pop(all_cell_data, laser_centers, sacc_centers, eye_grid, ...
        laser_bl_mask, sacc_bl_mask, eye_bl_mask, nBins_laser, nBins_sacc, nBins_eye);
    do_plot(all_cell_data, pop_data, laser_centers, sacc_centers, eye_grid, ...
        laser_window, saccade_window, ...
        color_vor_pref, color_vor_nonpref, color_spont_pref, color_spont_nonpref, ...
        shade_vor_pref, shade_vor_nonpref, shade_spont_pref, shade_spont_nonpref, ...
        color_laser_line, color_laser_bg, color_eye_line, color_eye_shade, shade_alpha, ...
        FontSize, TitleFontSize, TickFontSize, LineWidth_FR, ...
        doShow, doSavePNG, doSaveFIG, doSaveEPS, out_folder);
    if doSavePNG || doSaveFIG || doSaveEPS
        save(fullfile(out_folder, 'optotag_pop_data.mat'), 'all_cell_data', 'pop_data', 'params');
    end
    return
end

%% =====================================================================
%  FULL RUN
%  =====================================================================
all_cell_data = cell(nCells, 1);
fprintf('\n  Optotagged Population PETHs v3.3 — %d cells\n', nCells);

for iC = 1:nCells
    neuron_name = opto_neurons{iC};
    SSN = opto_sessions{iC};
    mouse_id = SSN(1:4);
    sess_path = fullfile(data_root, mouse_id, SSN);
    fprintf('\n--- %d/%d: %s ---\n', iC, nCells, neuron_name);

    cd_s = struct();
    cd_s.neuron_name = neuron_name;
    cd_s.SSN = SSN;

    try
        pushdir(sess_path);

        % Spikes
        [S_all, spike_cell_names] = LoadSpikesQuick();
        cell_idx = find_cell_index(spike_cell_names, neuron_name, SSN);
        if isempty(cell_idx)
            cd_s.status = 'ERROR: cell not found';
            all_cell_data{iC} = cd_s; popdir; continue
        end
        myCell = SelectTS([], S_all, cell_idx);

        % Laser
        [~, ~, laser_on, ~, n_pulses] = SortBrainstemEventLabels3;
        if isnan(laser_on)
            cd_s.status = 'ERROR: no laser';
            all_cell_data{iC} = cd_s; popdir; continue
        end
        cd_s.n_laser_pulses = n_pulses;
        cd_s.laser_dur = 1.0;
        if contains(SSN, 'M271'); cd_s.laser_dur = 0.5; end

        % Laser FR PETH
        [cd_s.laser_fr, cd_s.laser_fr_sem] = compute_peth_smooth(myCell, laser_on, laser_bins, ...
            peth_dt, peth_smooth_win, peth_excessBounds);

        % --- Laser eye position PETH ---
        sacc_file = strcat(SSN, '-saccades-edited.mat');
        has_eye = false;
        if exist(sacc_file, 'file')
            eye_vars = load(sacc_file, 'tsdH', 'k');
            if isfield(eye_vars, 'tsdH')
                has_eye = true;
                ep_tvec = eye_vars.tsdH.tvec(:);
                ep_data = eye_vars.tsdH.data(:);
                % Convert pixels to degrees using calibration factor k
                if isfield(eye_vars, 'k')
                    cal_k = eye_vars.k;
                else
                    cal_k = 1;
                    fprintf('  WARNING: no calibration factor k, eye pos in pixels\n');
                end
                ep_data = ep_data * cal_k;
            end
        end

        if has_eye
            eye_traces = nan(n_pulses, nBins_eye);
            for iP = 1:n_pulses
                t0 = laser_on(iP);
                search_mask = ep_tvec >= (t0 + laser_window(1) - 0.1) & ...
                              ep_tvec <= (t0 + laser_window(2) + 0.1);
                if sum(search_mask) < 5; continue; end
                rel_times = ep_tvec(search_mask) - t0;
                eye_traces(iP, :) = interp1(rel_times, ep_data(search_mask), eye_grid, 'linear', NaN);
            end
            % Remove all-NaN trials
            valid_trials = ~all(isnan(eye_traces), 2);
            eye_traces = eye_traces(valid_trials, :);

            % Baseline-subtract each trial
            for iT = 1:size(eye_traces, 1)
                bl_val = mean(eye_traces(iT, eye_bl_mask), 'omitnan');
                eye_traces(iT, :) = eye_traces(iT, :) - bl_val;
            end

            cd_s.eye_pos_mean = mean(eye_traces, 1, 'omitnan');
            cd_s.eye_pos_sem  = std(eye_traces, 0, 1, 'omitnan') / sqrt(size(eye_traces, 1));
            cd_s.eye_pos_traces = eye_traces;
            cd_s.n_eye_trials = size(eye_traces, 1);
            cd_s.has_eye = true;
            fprintf('  Eye position: %d valid trials\n', cd_s.n_eye_trials);
        else
            cd_s.eye_pos_mean = nan(1, nBins_eye);
            cd_s.eye_pos_sem  = nan(1, nBins_eye);
            cd_s.has_eye = false;
            fprintf('  No eye position data\n');
        end

        % --- Saccade times (Chenhao) + VOR/stationary split ---
        [temporal_all, nasal_all] = load_chenhao_saccades(SSN);

        stat_file = strcat(SSN, '-AHV_StationaryTimes.mat');
        if exist(stat_file, 'file')
            st = load(stat_file, 'STtstart', 'STtend');
            STtstart = st.STtstart(:); STtend = st.STtend(:);
            t_stat_mask = false(size(temporal_all));
            n_stat_mask = false(size(nasal_all));
            for iST = 1:length(STtstart)
                t_stat_mask = t_stat_mask | (temporal_all > STtstart(iST) & temporal_all < STtend(iST));
                n_stat_mask = n_stat_mask | (nasal_all > STtstart(iST) & nasal_all < STtend(iST));
            end
            temporal_move = temporal_all(~t_stat_mask);
            nasal_move    = nasal_all(~n_stat_mask);
            temporal_stat = temporal_all(t_stat_mask);
            nasal_stat    = nasal_all(n_stat_mask);
        else
            temporal_move = temporal_all; nasal_move = nasal_all;
            temporal_stat = []; nasal_stat = [];
        end

        % VOR saccade PETHs
        [cd_s.vor_fr_CCW, cd_s.vor_sem_CCW] = compute_peth_smooth(myCell, temporal_move, sacc_bins, peth_dt, peth_smooth_win, peth_excessBounds);
        [cd_s.vor_fr_CW, cd_s.vor_sem_CW]   = compute_peth_smooth(myCell, nasal_move, sacc_bins, peth_dt, peth_smooth_win, peth_excessBounds);
        cd_s.n_vor_CCW = length(temporal_move); cd_s.n_vor_CW = length(nasal_move);

        % Spontaneous saccade PETHs
        [cd_s.stat_fr_CCW, cd_s.stat_sem_CCW] = compute_peth_smooth(myCell, temporal_stat, sacc_bins, peth_dt, peth_smooth_win, peth_excessBounds);
        [cd_s.stat_fr_CW, cd_s.stat_sem_CW]   = compute_peth_smooth(myCell, nasal_stat, sacc_bins, peth_dt, peth_smooth_win, peth_excessBounds);
        cd_s.n_stat_CCW = length(temporal_stat); cd_s.n_stat_CW = length(nasal_stat);

        % Preferred direction (from VOR burst modulation)
        bl_CCW = mean(cd_s.vor_fr_CCW(sacc_bl_mask), 'omitnan');
        bst_CCW = mean(cd_s.vor_fr_CCW(sacc_burst_mask), 'omitnan');
        bl_CW = mean(cd_s.vor_fr_CW(sacc_bl_mask), 'omitnan');
        bst_CW = mean(cd_s.vor_fr_CW(sacc_burst_mask), 'omitnan');
        if bl_CCW > 0; bm_CCW = (bst_CCW - bl_CCW)/bl_CCW; else; bm_CCW = NaN; end
        if bl_CW > 0;  bm_CW  = (bst_CW - bl_CW)/bl_CW;   else; bm_CW = NaN;  end
        if isnan(bm_CW) || bm_CCW >= bm_CW; pref_dir = 'CCW'; else; pref_dir = 'CW'; end
        cd_s.pref_direction = pref_dir;

        if strcmp(pref_dir, 'CCW')
            cd_s.vor_pref = cd_s.vor_fr_CCW;   cd_s.vor_nonpref = cd_s.vor_fr_CW;
            cd_s.stat_pref = cd_s.stat_fr_CCW;  cd_s.stat_nonpref = cd_s.stat_fr_CW;
            cd_s.vor_pref_sem = cd_s.vor_sem_CCW;   cd_s.vor_nonpref_sem = cd_s.vor_sem_CW;
            cd_s.stat_pref_sem = cd_s.stat_sem_CCW;  cd_s.stat_nonpref_sem = cd_s.stat_sem_CW;
            cd_s.n_vor_pref = cd_s.n_vor_CCW;   cd_s.n_vor_nonpref = cd_s.n_vor_CW;
            cd_s.n_stat_pref = cd_s.n_stat_CCW;  cd_s.n_stat_nonpref = cd_s.n_stat_CW;
        else
            cd_s.vor_pref = cd_s.vor_fr_CW;     cd_s.vor_nonpref = cd_s.vor_fr_CCW;
            cd_s.stat_pref = cd_s.stat_fr_CW;    cd_s.stat_nonpref = cd_s.stat_fr_CCW;
            cd_s.vor_pref_sem = cd_s.vor_sem_CW;     cd_s.vor_nonpref_sem = cd_s.vor_sem_CCW;
            cd_s.stat_pref_sem = cd_s.stat_sem_CW;    cd_s.stat_nonpref_sem = cd_s.stat_sem_CCW;
            cd_s.n_vor_pref = cd_s.n_vor_CW;     cd_s.n_vor_nonpref = cd_s.n_vor_CCW;
            cd_s.n_stat_pref = cd_s.n_stat_CW;    cd_s.n_stat_nonpref = cd_s.n_stat_CCW;
        end

        fprintf('  Pref:%s | VOR: CCW=%d CW=%d | Stat: CCW=%d CW=%d | Laser: %d\n', ...
            pref_dir, cd_s.n_vor_CCW, cd_s.n_vor_CW, cd_s.n_stat_CCW, cd_s.n_stat_CW, n_pulses);

        cd_s.status = 'OK';
        popdir;
    catch ME
        cd_s.status = sprintf('ERROR: %s', ME.message);
        fprintf('  ERROR: %s\n', ME.message);
        try popdir; catch; end
    end
    all_cell_data{iC} = cd_s;
end

%% Population + plot
pop_data = compute_pop(all_cell_data, laser_centers, sacc_centers, eye_grid, ...
    laser_bl_mask, sacc_bl_mask, eye_bl_mask, nBins_laser, nBins_sacc, nBins_eye);

do_plot(all_cell_data, pop_data, laser_centers, sacc_centers, eye_grid, ...
    laser_window, saccade_window, ...
    color_vor_pref, color_vor_nonpref, color_spont_pref, color_spont_nonpref, ...
    shade_vor_pref, shade_vor_nonpref, shade_spont_pref, shade_spont_nonpref, ...
    color_laser_line, color_laser_bg, color_eye_line, color_eye_shade, shade_alpha, ...
    FontSize, TitleFontSize, TickFontSize, LineWidth_FR, ...
    doShow, doSavePNG, doSaveFIG, doSaveEPS, out_folder);

if doSavePNG || doSaveFIG || doSaveEPS
    save(fullfile(out_folder, 'optotag_pop_data.mat'), 'all_cell_data', 'pop_data', 'params');
    fprintf('\nSaved to: %s\n', out_folder);
end
end

%% =====================================================================
%  COMPUTE POPULATION DATA
%  =====================================================================
function pop = compute_pop(all_cell_data, laser_centers, sacc_centers, eye_grid, ...
        laser_bl_mask, sacc_bl_mask, eye_bl_mask, nBins_laser, nBins_sacc, nBins_eye)

    nCells = length(all_cell_data);
    vor_pref_mat = nan(nCells, nBins_sacc); vor_nonpref_mat = nan(nCells, nBins_sacc);
    stat_pref_mat = nan(nCells, nBins_sacc); stat_nonpref_mat = nan(nCells, nBins_sacc);
    laser_mat = nan(nCells, nBins_laser);
    eye_pos_mat = nan(nCells, nBins_eye);
    valid = false(nCells, 1);
    eye_valid = false(nCells, 1);

    for iC = 1:nCells
        cd = all_cell_data{iC};
        if ~isstruct(cd) || ~strcmp(cd.status, 'OK'); continue; end
        valid(iC) = true;
        vor_pref_mat(iC,:) = cd.vor_pref;
        vor_nonpref_mat(iC,:) = cd.vor_nonpref;
        stat_pref_mat(iC,:) = cd.stat_pref;
        stat_nonpref_mat(iC,:) = cd.stat_nonpref;
        laser_mat(iC,:) = cd.laser_fr;
        if isfield(cd, 'has_eye') && cd.has_eye
            eye_pos_mat(iC,:) = cd.eye_pos_mean;
            eye_valid(iC) = true;
        end
    end

    n_valid = sum(valid);
    n_eye_valid = sum(eye_valid);
    pop = struct();
    pop.n_valid = n_valid;
    pop.n_eye_valid = n_eye_valid;

    % Raw mean/SEM
    pop.raw_vor_pref_mean     = mean(vor_pref_mat(valid,:), 1, 'omitnan');
    pop.raw_vor_pref_sem      = std(vor_pref_mat(valid,:), 0, 1, 'omitnan') / sqrt(n_valid);
    pop.raw_vor_nonpref_mean  = mean(vor_nonpref_mat(valid,:), 1, 'omitnan');
    pop.raw_vor_nonpref_sem   = std(vor_nonpref_mat(valid,:), 0, 1, 'omitnan') / sqrt(n_valid);
    pop.raw_stat_pref_mean    = mean(stat_pref_mat(valid,:), 1, 'omitnan');
    pop.raw_stat_pref_sem     = std(stat_pref_mat(valid,:), 0, 1, 'omitnan') / sqrt(n_valid);
    pop.raw_stat_nonpref_mean = mean(stat_nonpref_mat(valid,:), 1, 'omitnan');
    pop.raw_stat_nonpref_sem  = std(stat_nonpref_mat(valid,:), 0, 1, 'omitnan') / sqrt(n_valid);
    pop.raw_laser_mean        = mean(laser_mat(valid,:), 1, 'omitnan');
    pop.raw_laser_sem         = std(laser_mat(valid,:), 0, 1, 'omitnan') / sqrt(n_valid);

    % Eye position population (already baseline-subtracted per trial)
    if n_eye_valid >= 2
        pop.eye_pos_mean = mean(eye_pos_mat(eye_valid,:), 1, 'omitnan');
        pop.eye_pos_sem  = std(eye_pos_mat(eye_valid,:), 0, 1, 'omitnan') / sqrt(n_eye_valid);
    else
        pop.eye_pos_mean = nan(1, nBins_eye);
        pop.eye_pos_sem  = nan(1, nBins_eye);
    end

    % Z-scored
    z_vor_pref = nan(nCells, nBins_sacc); z_vor_nonpref = nan(nCells, nBins_sacc);
    z_stat_pref = nan(nCells, nBins_sacc); z_stat_nonpref = nan(nCells, nBins_sacc);
    z_laser = nan(nCells, nBins_laser);
    for iC = 1:nCells
        if ~valid(iC); continue; end
        s_bl = mean(vor_pref_mat(iC, sacc_bl_mask), 'omitnan');
        s_sd = std(vor_pref_mat(iC, sacc_bl_mask), 'omitnan');
        if s_sd < 1; s_sd = 1; end
        z_vor_pref(iC,:)     = (vor_pref_mat(iC,:) - s_bl) / s_sd;
        z_vor_nonpref(iC,:)  = (vor_nonpref_mat(iC,:) - s_bl) / s_sd;
        z_stat_pref(iC,:)    = (stat_pref_mat(iC,:) - s_bl) / s_sd;
        z_stat_nonpref(iC,:) = (stat_nonpref_mat(iC,:) - s_bl) / s_sd;
        l_bl = mean(laser_mat(iC, laser_bl_mask), 'omitnan');
        l_sd = std(laser_mat(iC, laser_bl_mask), 'omitnan');
        if l_sd < 1; l_sd = 1; end
        z_laser(iC,:) = (laser_mat(iC,:) - l_bl) / l_sd;
    end
    pop.z_vor_pref_mean     = mean(z_vor_pref(valid,:), 1, 'omitnan');
    pop.z_vor_pref_sem      = std(z_vor_pref(valid,:), 0, 1, 'omitnan') / sqrt(n_valid);
    pop.z_vor_nonpref_mean  = mean(z_vor_nonpref(valid,:), 1, 'omitnan');
    pop.z_vor_nonpref_sem   = std(z_vor_nonpref(valid,:), 0, 1, 'omitnan') / sqrt(n_valid);
    pop.z_stat_pref_mean    = mean(z_stat_pref(valid,:), 1, 'omitnan');
    pop.z_stat_pref_sem     = std(z_stat_pref(valid,:), 0, 1, 'omitnan') / sqrt(n_valid);
    pop.z_stat_nonpref_mean = mean(z_stat_nonpref(valid,:), 1, 'omitnan');
    pop.z_stat_nonpref_sem  = std(z_stat_nonpref(valid,:), 0, 1, 'omitnan') / sqrt(n_valid);
    pop.z_laser_mean        = mean(z_laser(valid,:), 1, 'omitnan');
    pop.z_laser_sem         = std(z_laser(valid,:), 0, 1, 'omitnan') / sqrt(n_valid);
end

%% =====================================================================
%  PLOT ALL
%  =====================================================================
function do_plot(all_cell_data, pop_data, laser_centers, sacc_centers, eye_grid, ...
        laser_window, saccade_window, ...
        c_vp, c_vn, c_sp, c_sn, sh_vp, sh_vn, sh_sp, sh_sn, ...
        c_laser, c_laser_bg, c_eye, c_eye_sh, sh_alpha, ...
        FS, TFS, TkFS, LW, doShow, doSavePNG, doSaveFIG, doSaveEPS, out_folder)

    %% Per-neuron 2x3
    for iC = 1:length(all_cell_data)
        cd = all_cell_data{iC};
        if ~isstruct(cd) || ~strcmp(cd.status, 'OK'); continue; end

        if doShow; fh = figure('Position', [50 50 1500 700]);
        else; fh = figure('Position', [50 50 1500 700], 'Visible', 'off'); end
        fh.UserData.date_created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        fh.UserData.generating_function = 'plot_optotag_population_PETHs v3.3';
        laser_dur = cd.laser_dur;

        % Helper for SEM fill
        sem_f = @(ax_h, x_v, m_v, s_v, c_sh) fill(ax_h, [x_v fliplr(x_v)], ...
            [m_v+s_v fliplr(m_v-s_v)], c_sh, 'EdgeColor', 'none', ...
            'FaceAlpha', sh_alpha, 'HandleVisibility', 'off');

        % (1,1) VOR Preferred
        ax1 = subplot(2,3,1); hold(ax1,'on');
        if ~any(isnan(cd.vor_pref + cd.vor_pref_sem))
            sem_f(ax1, sacc_centers, cd.vor_pref, cd.vor_pref_sem, sh_vp);
        end
        plot(ax1, sacc_centers, cd.vor_pref, '-', 'Color', c_vp, 'LineWidth', LW);
        xline(ax1, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1.5);
        ylabel(ax1, 'FR (Hz)', 'FontSize', FS);
        title(ax1, sprintf('VOR — Preferred (n=%d)', cd.n_vor_pref), 'FontSize', TFS, 'FontWeight', 'bold');
        xlim(ax1, saccade_window); set(ax1, 'FontSize', TkFS); box(ax1,'off');

        % (1,2) VOR Non-Preferred
        ax2 = subplot(2,3,2); hold(ax2,'on');
        if ~any(isnan(cd.vor_nonpref + cd.vor_nonpref_sem))
            sem_f(ax2, sacc_centers, cd.vor_nonpref, cd.vor_nonpref_sem, sh_vn);
        end
        plot(ax2, sacc_centers, cd.vor_nonpref, '-', 'Color', c_vn, 'LineWidth', LW);
        xline(ax2, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1.5);
        title(ax2, sprintf('VOR — Non-Preferred (n=%d)', cd.n_vor_nonpref), 'FontSize', TFS, 'FontWeight', 'bold');
        xlim(ax2, saccade_window); set(ax2, 'FontSize', TkFS); box(ax2,'off');

        % (1,3) Laser FR
        ax3 = subplot(2,3,3); hold(ax3,'on');
        yl3_max = max(cd.laser_fr + cd.laser_fr_sem);
        if isnan(yl3_max) || yl3_max == 0; yl3_max = 10; end
        yl3 = [0 yl3_max*1.15];
        fill(ax3, [0 laser_dur laser_dur 0], [yl3(1) yl3(1) yl3(2) yl3(2)], ...
            c_laser_bg, 'EdgeColor', 'none', 'FaceAlpha', 0.25, 'HandleVisibility', 'off');
        if ~any(isnan(cd.laser_fr + cd.laser_fr_sem))
            fill(ax3, [laser_centers fliplr(laser_centers)], ...
                [cd.laser_fr+cd.laser_fr_sem fliplr(cd.laser_fr-cd.laser_fr_sem)], ...
                [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
        end
        plot(ax3, laser_centers, cd.laser_fr, '-', 'Color', c_laser, 'LineWidth', LW);
        xline(ax3, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1.5);
        title(ax3, sprintf('Laser FR (%d pulses)', cd.n_laser_pulses), 'FontSize', TFS, 'FontWeight', 'bold');
        xlim(ax3, laser_window); ylim(ax3, yl3); set(ax3, 'FontSize', TkFS); box(ax3,'off');

        % (2,1) Spontaneous Preferred
        ax4 = subplot(2,3,4); hold(ax4,'on');
        if ~any(isnan(cd.stat_pref + cd.stat_pref_sem))
            sem_f(ax4, sacc_centers, cd.stat_pref, cd.stat_pref_sem, sh_sp);
        end
        plot(ax4, sacc_centers, cd.stat_pref, '-', 'Color', c_sp, 'LineWidth', LW);
        xline(ax4, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1.5);
        xlabel(ax4, 'Time (s)', 'FontSize', FS); ylabel(ax4, 'FR (Hz)', 'FontSize', FS);
        title(ax4, sprintf('Spont — Preferred (n=%d)', cd.n_stat_pref), 'FontSize', TFS, 'FontWeight', 'bold');
        xlim(ax4, saccade_window); set(ax4, 'FontSize', TkFS); box(ax4,'off');

        % (2,2) Spontaneous Non-Preferred
        ax5 = subplot(2,3,5); hold(ax5,'on');
        if ~any(isnan(cd.stat_nonpref + cd.stat_nonpref_sem))
            sem_f(ax5, sacc_centers, cd.stat_nonpref, cd.stat_nonpref_sem, sh_sn);
        end
        plot(ax5, sacc_centers, cd.stat_nonpref, '-', 'Color', c_sn, 'LineWidth', LW);
        xline(ax5, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1.5);
        xlabel(ax5, 'Time (s)', 'FontSize', FS);
        title(ax5, sprintf('Spont — Non-Preferred (n=%d)', cd.n_stat_nonpref), 'FontSize', TFS, 'FontWeight', 'bold');
        xlim(ax5, saccade_window); set(ax5, 'FontSize', TkFS); box(ax5,'off');

        % (2,3) Laser Eye Position
        ax6 = subplot(2,3,6); hold(ax6,'on');
        if isfield(cd, 'has_eye') && cd.has_eye
            % Blue background for laser-on
            yl6_data = [cd.eye_pos_mean - cd.eye_pos_sem, cd.eye_pos_mean + cd.eye_pos_sem];
            yl6 = [min(yl6_data)*1.2 max(yl6_data)*1.2];
            if diff(yl6) < 0.5; yl6 = [-5 5]; end
            fill(ax6, [0 laser_dur laser_dur 0], [yl6(1) yl6(1) yl6(2) yl6(2)], ...
                c_laser_bg, 'EdgeColor', 'none', 'FaceAlpha', 0.25, 'HandleVisibility', 'off');
            % Individual traces (light)
            if isfield(cd, 'eye_pos_traces')
                for iT = 1:min(size(cd.eye_pos_traces, 1), 50)
                    plot(ax6, eye_grid, cd.eye_pos_traces(iT,:), '-', ...
                        'Color', [0.75 0.75 0.75 0.3], 'LineWidth', 0.3);
                end
            end
            % Mean
            plot(ax6, eye_grid, cd.eye_pos_mean, '-', 'Color', c_eye, 'LineWidth', LW);
            xline(ax6, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1.5);
            ylim(ax6, yl6);
        else
            text(ax6, 0.5, 0.5, 'No eye data', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', FS, 'Color', [.6 .6 .6]);
        end
        xlabel(ax6, 'Time (s)', 'FontSize', FS);
        ylabel(ax6, 'Eye position (deg)', 'FontSize', FS);
        title(ax6, sprintf('Laser — Eye Position (%d pulses)', cd.n_laser_pulses), 'FontSize', TFS, 'FontWeight', 'bold');
        xlim(ax6, laser_window); set(ax6, 'FontSize', TkFS); box(ax6,'off');

        % Match saccade y-axes per row
        yl1 = ylim(ax1); yl2 = ylim(ax2);
        yl_r1 = [min(yl1(1), yl2(1)) max(yl1(2), yl2(2))];
        ylim(ax1, yl_r1); ylim(ax2, yl_r1);
        yl4 = ylim(ax4); yl5 = ylim(ax5);
        yl_r2 = [min(yl4(1), yl5(1)) max(yl4(2), yl5(2))];
        ylim(ax4, yl_r2); ylim(ax5, yl_r2);

        sgtitle(sprintf('%s (Pref: %s)', cd.neuron_name, cd.pref_direction), ...
            'FontSize', TFS+4, 'FontWeight', 'bold', 'Interpreter', 'none');

        safe_name = strrep(cd.neuron_name, '-', '_');
        save_fig_trio(fh, fullfile(out_folder, sprintf('optotag_%s', safe_name)), ...
            doSavePNG, doSaveFIG, doSaveEPS, doShow);
    end

    %% Population raw FR (2x3)
    plot_pop_2x3(pop_data, 'raw', laser_centers, sacc_centers, eye_grid, ...
        laser_window, saccade_window, ...
        c_vp, c_vn, c_sp, c_sn, sh_vp, sh_vn, sh_sp, sh_sn, ...
        c_laser, c_laser_bg, c_eye, c_eye_sh, sh_alpha, ...
        FS, TFS, TkFS, LW, doShow, doSavePNG, doSaveFIG, doSaveEPS, out_folder);

    %% Population z-scored (2x3)
    plot_pop_2x3(pop_data, 'zscore', laser_centers, sacc_centers, eye_grid, ...
        laser_window, saccade_window, ...
        c_vp, c_vn, c_sp, c_sn, sh_vp, sh_vn, sh_sp, sh_sn, ...
        c_laser, c_laser_bg, c_eye, c_eye_sh, sh_alpha, ...
        FS, TFS, TkFS, LW, doShow, doSavePNG, doSaveFIG, doSaveEPS, out_folder);
end

%% =====================================================================
%  POPULATION 2x3 FIGURE
%  =====================================================================
function plot_pop_2x3(pop, mode, lc, sc, eg, lw_ax, sw_ax, ...
        c_vp, c_vn, c_sp, c_sn, sh_vp, sh_vn, sh_sp, sh_sn, ...
        c_laser, c_laser_bg, c_eye, c_eye_sh, sh_alpha, ...
        FS, TFS, TkFS, LW, doShow, doSavePNG, doSaveFIG, doSaveEPS, out_folder)

    if strcmp(mode, 'zscore')
        vp_m = pop.z_vor_pref_mean; vp_s = pop.z_vor_pref_sem;
        vn_m = pop.z_vor_nonpref_mean; vn_s = pop.z_vor_nonpref_sem;
        sp_m = pop.z_stat_pref_mean; sp_s = pop.z_stat_pref_sem;
        sn_m = pop.z_stat_nonpref_mean; sn_s = pop.z_stat_nonpref_sem;
        lr_m = pop.z_laser_mean; lr_s = pop.z_laser_sem;
        y_label = 'FR (z-score)'; fname = 'population_optotag_zscore';
    else
        vp_m = pop.raw_vor_pref_mean; vp_s = pop.raw_vor_pref_sem;
        vn_m = pop.raw_vor_nonpref_mean; vn_s = pop.raw_vor_nonpref_sem;
        sp_m = pop.raw_stat_pref_mean; sp_s = pop.raw_stat_pref_sem;
        sn_m = pop.raw_stat_nonpref_mean; sn_s = pop.raw_stat_nonpref_sem;
        lr_m = pop.raw_laser_mean; lr_s = pop.raw_laser_sem;
        y_label = 'FR (Hz)'; fname = 'population_optotag_raw';
    end

    if doShow; fh = figure('Position', [50 50 1500 700]);
    else; fh = figure('Position', [50 50 1500 700], 'Visible', 'off'); end
    fh.UserData.date_created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fh.UserData.generating_function = 'plot_optotag_population_PETHs v3.3';

    sem_fill = @(ax, x, m, s, c_sh) fill(ax, [x fliplr(x)], [m+s fliplr(m-s)], ...
        c_sh, 'EdgeColor', 'none', 'FaceAlpha', sh_alpha, 'HandleVisibility', 'off');

    % (1,1) VOR Preferred
    ax1 = subplot(2,3,1); hold(ax1,'on');
    if ~any(isnan(vp_m+vp_s)); sem_fill(ax1, sc, vp_m, vp_s, sh_vp); end
    plot(ax1, sc, vp_m, '-', 'Color', c_vp, 'LineWidth', LW);
    xline(ax1, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1.5);
    ylabel(ax1, y_label, 'FontSize', FS);
    title(ax1, 'VOR — Preferred', 'FontSize', TFS, 'FontWeight', 'bold');
    xlim(ax1, sw_ax); set(ax1, 'FontSize', TkFS); box(ax1,'off');

    % (1,2) VOR Non-Preferred
    ax2 = subplot(2,3,2); hold(ax2,'on');
    if ~any(isnan(vn_m+vn_s)); sem_fill(ax2, sc, vn_m, vn_s, sh_vn); end
    plot(ax2, sc, vn_m, '-', 'Color', c_vn, 'LineWidth', LW);
    xline(ax2, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1.5);
    title(ax2, 'VOR — Non-Preferred', 'FontSize', TFS, 'FontWeight', 'bold');
    xlim(ax2, sw_ax); set(ax2, 'FontSize', TkFS); box(ax2,'off');

    % (1,3) Laser FR
    ax3 = subplot(2,3,3); hold(ax3,'on');
    all_lr = [lr_m - lr_s, lr_m + lr_s];
    yl_lr = [min(all_lr)*0.9 max(all_lr)*1.1]; if diff(yl_lr) == 0; yl_lr = [0 10]; end
    fill(ax3, [0 1 1 0], [yl_lr(1) yl_lr(1) yl_lr(2) yl_lr(2)], ...
        c_laser_bg, 'EdgeColor', 'none', 'FaceAlpha', 0.25, 'HandleVisibility', 'off');
    if ~any(isnan(lr_m+lr_s))
        fill(ax3, [lc fliplr(lc)], [lr_m+lr_s fliplr(lr_m-lr_s)], ...
            [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
    end
    plot(ax3, lc, lr_m, '-', 'Color', c_laser, 'LineWidth', LW);
    xline(ax3, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1.5);
    title(ax3, 'Laser — FR', 'FontSize', TFS, 'FontWeight', 'bold');
    xlim(ax3, lw_ax); ylim(ax3, yl_lr); set(ax3, 'FontSize', TkFS); box(ax3,'off');

    % (2,1) Spontaneous Preferred
    ax4 = subplot(2,3,4); hold(ax4,'on');
    if ~any(isnan(sp_m+sp_s)); sem_fill(ax4, sc, sp_m, sp_s, sh_sp); end
    plot(ax4, sc, sp_m, '-', 'Color', c_sp, 'LineWidth', LW);
    xline(ax4, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1.5);
    xlabel(ax4, 'Time (s)', 'FontSize', FS); ylabel(ax4, y_label, 'FontSize', FS);
    title(ax4, 'Spont — Preferred', 'FontSize', TFS, 'FontWeight', 'bold');
    xlim(ax4, sw_ax); set(ax4, 'FontSize', TkFS); box(ax4,'off');

    % (2,2) Spontaneous Non-Preferred
    ax5 = subplot(2,3,5); hold(ax5,'on');
    if ~any(isnan(sn_m+sn_s)); sem_fill(ax5, sc, sn_m, sn_s, sh_sn); end
    plot(ax5, sc, sn_m, '-', 'Color', c_sn, 'LineWidth', LW);
    xline(ax5, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1.5);
    xlabel(ax5, 'Time (s)', 'FontSize', FS);
    title(ax5, 'Spont — Non-Preferred', 'FontSize', TFS, 'FontWeight', 'bold');
    xlim(ax5, sw_ax); set(ax5, 'FontSize', TkFS); box(ax5,'off');

    % (2,3) Laser Eye Position
    ax6 = subplot(2,3,6); hold(ax6,'on');
    ep_m = pop.eye_pos_mean; ep_s = pop.eye_pos_sem;
    if ~all(isnan(ep_m))
        all_ep = [ep_m - ep_s, ep_m + ep_s];
        yl_ep = [min(all_ep)*1.2 max(all_ep)*1.2]; if diff(yl_ep) < 0.5; yl_ep = [-5 5]; end
        fill(ax6, [0 1 1 0], [yl_ep(1) yl_ep(1) yl_ep(2) yl_ep(2)], ...
            c_laser_bg, 'EdgeColor', 'none', 'FaceAlpha', 0.25, 'HandleVisibility', 'off');
        if ~any(isnan(ep_m + ep_s))
            sem_fill(ax6, eg, ep_m, ep_s, c_eye_sh);
        end
        plot(ax6, eg, ep_m, '-', 'Color', c_eye, 'LineWidth', LW);
        xline(ax6, 0, '--', 'Color', [.3 .3 .3], 'LineWidth', 1.5);
        ylim(ax6, yl_ep);
    end
    xlabel(ax6, 'Time (s)', 'FontSize', FS);
    ylabel(ax6, 'Eye position (deg)', 'FontSize', FS);
    title(ax6, sprintf('Laser — Eye Position (%d cells)', pop.n_eye_valid), 'FontSize', TFS, 'FontWeight', 'bold');
    xlim(ax6, lw_ax); set(ax6, 'FontSize', TkFS); box(ax6,'off');

    % Match saccade y-axes per row
    yl1 = ylim(ax1); yl2 = ylim(ax2);
    yl_r1 = [min(yl1(1), yl2(1)) max(yl1(2), yl2(2))];
    ylim(ax1, yl_r1); ylim(ax2, yl_r1);
    yl4 = ylim(ax4); yl5 = ylim(ax5);
    yl_r2 = [min(yl4(1), yl5(1)) max(yl4(2), yl5(2))];
    ylim(ax4, yl_r2); ylim(ax5, yl_r2);

    if strcmp(mode, 'zscore'); ns = 'z-scored'; else; ns = 'raw FR'; end
    sgtitle(sprintf('Optotagged Population — %s (%d cells)', ns, pop.n_valid), ...
        'FontSize', TFS+4, 'FontWeight', 'bold');
    save_fig_trio(fh, fullfile(out_folder, fname), doSavePNG, doSaveFIG, doSaveEPS, doShow);
end

%% =====================================================================
%  HELPERS
%  =====================================================================
function [fr_smooth, fr_sem] = compute_peth_smooth(myCell, event_times, peth_bins, peth_dt, smooth_win, excessBounds)
    nBins = length(peth_bins) - 1;
    nEvents = length(event_times);
    if nEvents < 1
        fr_smooth = nan(1, nBins); fr_sem = nan(1, nBins); return
    end
    pw = [peth_bins(1) peth_bins(end)];

    % Per-trial FR matrix for SEM
    trial_fr = zeros(nEvents, nBins);
    spk_rel = [];
    for iE = 1:nEvents
        S0 = restrict(myCell, event_times(iE)+pw(1)-excessBounds, event_times(iE)+pw(2)+excessBounds);
        S0 = restrict(S0, event_times(iE)+pw(1), event_times(iE)+pw(2));
        if ~isempty(S0.t{1})
            spk_rel = [spk_rel; S0.t{1} - event_times(iE)]; %#ok<AGROW>
            trial_hist = histcounts(S0.t{1} - event_times(iE), peth_bins);
            trial_fr(iE, :) = trial_hist / peth_dt;
        end
    end

    % Population-averaged PETH (same as before)
    m_hist = histcounts(spk_rel, peth_bins);
    fr_peth = m_hist / peth_dt / nEvents;
    fr_smooth = smoothdata(movmean(fr_peth, smooth_win), 'gaussian', smooth_win);

    % SEM across trials (smooth each trial first, then SEM)
    trial_fr_smooth = smoothdata(movmean(trial_fr, smooth_win, 2), 2, 'gaussian', smooth_win);
    fr_sem = std(trial_fr_smooth, 0, 1, 'omitnan') / sqrt(nEvents);
end

function cell_idx = find_cell_index(names, neuron, SSN)
    cell_idx = [];
    for cc = 1:length(names)
        if strcmp(names{cc}, neuron); cell_idx = cc; return; end
    end
    for cc = 1:length(names)
        if contains(names{cc}, strrep(neuron, [SSN '-'], '')); cell_idx = cc; return; end
    end
end

function save_fig_trio(fh, fp, doP, doF, doE, doS)
    if doP; print(fh, [fp '.png'], '-dpng', '-r300'); end
    if doF; savefig(fh, [fp '.fig']); end
    if doE; print(fh, [fp '.eps'], '-depsc2'); end
    if ~doS; close(fh); end
end
