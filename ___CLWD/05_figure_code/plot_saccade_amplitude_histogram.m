function [amp_data, params] = plot_saccade_amplitude_histogram(varargin)
% plot_saccade_amplitude_histogram  v1.0
% Histogram of saccade amplitudes: VOR (moving) vs Spontaneous (stationary).
% Pools across all sessions in NPH_all_tfilelist.
%
%   Uses pre-computed nasalAmplitudes/temporalAmplitudes from saccades-edited.mat
%   (multiplied by calibration factor k). Splits VOR vs stationary using
%   STtstart/STtend from AHV_StationaryTimes.mat. Uses Jeff's saccade times
%   for amplitude indexing (amplitudes are pre-computed per Jeff's events).
%
%   INPUTS (name-value pairs via varargin):
%     tfilelist_path    - path to NPH_all_tfilelist.mat
%     pixel_ceiling     - max amplitude in deg (default: 45)
%     bin_edges         - histogram bin edges in deg (default: 0:1:45)
%     doSavePNG/FIG/EPS - save flags (default: 1)
%     doShow            - display figure (default: 1)
%     save_dir          - output root
%
%   OUTPUTS:
%     amp_data - struct with all amplitude vectors and session info
%     params   - struct of parameter values
%
%   DEPENDENCIES:
%     HD_GetSSN, pushdir, popdir, process_varargin

%% =====================================================================
%  DEFAULTS
%  =====================================================================
tfilelist_path    = 'C:\Jeff\U01\datatouse\NPH all tfilelist.mat';
save_dir          = 'C:\Jeff\U01\datatouse\___analysis_output';
precomputed_amp   = [];   % pass amp_data struct from prior run to skip recomputation
pixel_ceiling     = 45;
amp_bin_edges     = 0:1:45;

% Figure
FontSize          = 20;
TitleFontSize     = 18;
TickFontSize      = 14;

% Output control
doSavePNG         = 1;
doSaveFIG         = 1;
doSaveEPS         = 1;
doShow            = 1;
verbose           = 1;

%% Process varargin
process_varargin(varargin);

% Manual extraction for struct params (process_varargin can silently fail on structs)
for iV = 1:2:length(varargin)
    if ischar(varargin{iV}) && strcmp(varargin{iV}, 'precomputed_amp')
        precomputed_amp = varargin{iV+1};
    end
end
%  =====================================================================
params = struct();
params.function_name    = 'plot_saccade_amplitude_histogram';
params.version          = '1.0';
params.date_run         = datestr(now, 'yyyy-mm-dd HH:MM:SS');
params.pixel_ceiling    = pixel_ceiling;
params.amp_bin_edges    = amp_bin_edges;
params.amplitude_source = 'nasalAmplitudes/temporalAmplitudes from saccades-edited.mat * k';
params.split_method     = 'STtstart/STtend from AHV_StationaryTimes.mat';

%% =====================================================================
%  WORKSPACE MODE
%  =====================================================================
if ~isempty(precomputed_amp) && isfield(precomputed_amp, 'all_VOR')
    fprintf('=== WORKSPACE MODE (amp_data from prior run) ===\n');
    amp_data = precomputed_amp;

    % Output folder
    run_date_str = datestr(now, 'yyyy-mm-dd');
    run_num = 1;
    out_folder = fullfile(save_dir, sprintf('saccade_amplitude_%s_run%d', run_date_str, run_num));
    while exist(out_folder, 'dir')
        run_num = run_num + 1;
        out_folder = fullfile(save_dir, sprintf('saccade_amplitude_%s_run%d', run_date_str, run_num));
    end
    if doSavePNG || doSaveFIG || doSaveEPS; mkdir(out_folder); end
    params.output_folder = out_folder;

    plot_amp_figure(amp_data, amp_bin_edges, FontSize, TitleFontSize, TickFontSize, ...
        doShow, doSavePNG, doSaveFIG, doSaveEPS, out_folder);

    if doSavePNG || doSaveFIG || doSaveEPS
        save(fullfile(out_folder, 'amplitude_data.mat'), 'amp_data', 'params');
        fprintf('Saved to: %s\n', out_folder);
    end
    return
end

%% =====================================================================
%  GET UNIQUE SESSIONS FROM TFILELIST
%  =====================================================================
tfl = load(tfilelist_path, 'NPH_all_tfilelist');
cell_list = tfl.NPH_all_tfilelist(:);

session_paths = cell(length(cell_list), 1);
for iN = 1:length(cell_list)
    session_paths{iN} = fileparts(cell_list{iN});
end
unique_sessions = unique(session_paths, 'stable');
nSessions = length(unique_sessions);
if verbose; fprintf('Found %d unique sessions from %d neurons.\n', nSessions, length(cell_list)); end

%% =====================================================================
%  OUTPUT FOLDER
%  =====================================================================
run_date_str = datestr(now, 'yyyy-mm-dd');
run_num = 1;
out_folder = fullfile(save_dir, sprintf('saccade_amplitude_%s_run%d', run_date_str, run_num));
while exist(out_folder, 'dir')
    run_num = run_num + 1;
    out_folder = fullfile(save_dir, sprintf('saccade_amplitude_%s_run%d', run_date_str, run_num));
end
if doSavePNG || doSaveFIG || doSaveEPS; mkdir(out_folder); end
params.output_folder = out_folder;

%% =====================================================================
%  COLLECT AMPLITUDES ACROSS ALL SESSIONS
%  =====================================================================
all_nasal_VOR     = [];
all_nasal_STAT    = [];
all_temporal_VOR  = [];
all_temporal_STAT = [];
sessions_used     = {};
sessions_skipped  = {};

for iS = 1:nSessions
    sess_path = unique_sessions{iS};

    try
        pushdir(sess_path);
        SSN = HD_GetSSN;

        % Load saccade data
        sacc_file = strcat(SSN, '-saccades-edited.mat');
        if ~exist(sacc_file, 'file')
            if verbose; fprintf('  %s: no saccades file, skipping\n', SSN); end
            sessions_skipped{end+1} = SSN; %#ok<AGROW>
            popdir; continue
        end

        sacc_data = load(sacc_file, 'nasalSaccades', 'temporalSaccades', ...
            'nasalAmplitudes', 'temporalAmplitudes', 'k');

        % Check for required fields
        if ~isfield(sacc_data, 'nasalAmplitudes') || ~isfield(sacc_data, 'temporalAmplitudes')
            if verbose; fprintf('  %s: no amplitude data, skipping\n', SSN); end
            sessions_skipped{end+1} = SSN; %#ok<AGROW>
            popdir; continue
        end

        % Get calibration factor
        if isfield(sacc_data, 'k')
            cal_k = sacc_data.k;
        else
            cal_k = 1;
            if verbose; fprintf('  %s: no calibration factor k, using k=1\n', SSN); end
        end

        % Get saccade times and amplitudes
        temp_times = sacc_data.temporalSaccades(:);
        nasal_times = sacc_data.nasalSaccades(:);
        temp_amps = sacc_data.temporalAmplitudes(:) * cal_k;
        nasal_amps = sacc_data.nasalAmplitudes(:) * cal_k;

        % Remove NaNs (maintaining index correspondence)
        temp_valid = ~isnan(temp_times) & ~isnan(temp_amps);
        nasal_valid = ~isnan(nasal_times) & ~isnan(nasal_amps);
        temp_times = temp_times(temp_valid);
        temp_amps  = temp_amps(temp_valid);
        nasal_times = nasal_times(nasal_valid);
        nasal_amps  = nasal_amps(nasal_valid);

        % Load stationary times
        stat_file = strcat(SSN, '-AHV_StationaryTimes.mat');
        if ~exist(stat_file, 'file')
            if verbose; fprintf('  %s: no StationaryTimes, skipping\n', SSN); end
            sessions_skipped{end+1} = SSN; %#ok<AGROW>
            popdir; continue
        end
        stat_info = load(stat_file, 'STtstart', 'STtend');
        STtstart = stat_info.STtstart(:);
        STtend   = stat_info.STtend(:);

        % Split temporal saccades: stationary vs VOR
        t_stat_mask = false(size(temp_times));
        for iST = 1:length(STtstart)
            t_stat_mask = t_stat_mask | (temp_times > STtstart(iST) & temp_times < STtend(iST));
        end

        % Split nasal saccades: stationary vs VOR
        n_stat_mask = false(size(nasal_times));
        for iST = 1:length(STtstart)
            n_stat_mask = n_stat_mask | (nasal_times > STtstart(iST) & nasal_times < STtend(iST));
        end

        % Apply pixel ceiling and collect
        t_vor  = abs(temp_amps(~t_stat_mask));
        t_stat = abs(temp_amps(t_stat_mask));
        n_vor  = abs(nasal_amps(~n_stat_mask));
        n_stat = abs(nasal_amps(n_stat_mask));

        t_vor  = t_vor(t_vor < pixel_ceiling);
        t_stat = t_stat(t_stat < pixel_ceiling);
        n_vor  = n_vor(n_vor < pixel_ceiling);
        n_stat = n_stat(n_stat < pixel_ceiling);

        all_temporal_VOR  = [all_temporal_VOR; t_vor]; %#ok<AGROW>
        all_temporal_STAT = [all_temporal_STAT; t_stat]; %#ok<AGROW>
        all_nasal_VOR     = [all_nasal_VOR; n_vor]; %#ok<AGROW>
        all_nasal_STAT    = [all_nasal_STAT; n_stat]; %#ok<AGROW>

        sessions_used{end+1} = SSN; %#ok<AGROW>
        popdir;

    catch ME
        if verbose; fprintf('  %s: ERROR %s\n', SSN, ME.message); end
        sessions_skipped{end+1} = SSN; %#ok<AGROW>
        try popdir; catch; end
    end
end

%% =====================================================================
%  COMBINE DIRECTIONS
%  =====================================================================
all_VOR  = [all_temporal_VOR; all_nasal_VOR];
all_STAT = [all_temporal_STAT; all_nasal_STAT];

fprintf('\n=== AMPLITUDE SUMMARY ===\n');
fprintf('  Sessions used:      %d\n', length(sessions_used));
fprintf('  Sessions skipped:   %d\n', length(sessions_skipped));
fprintf('  VOR saccades:       %d (temporal: %d, nasal: %d)\n', ...
    length(all_VOR), length(all_temporal_VOR), length(all_nasal_VOR));
fprintf('  Stationary saccades: %d (temporal: %d, nasal: %d)\n', ...
    length(all_STAT), length(all_temporal_STAT), length(all_nasal_STAT));
fprintf('  VOR amplitude:      mean=%.1f, median=%.1f deg\n', mean(all_VOR), median(all_VOR));
fprintf('  Stat amplitude:     mean=%.1f, median=%.1f deg\n', mean(all_STAT), median(all_STAT));
[p_rs, ~] = ranksum(all_VOR, all_STAT);
fprintf('  Ranksum p-value:    %.4e\n', p_rs);
fprintf('=========================\n');

%% =====================================================================
%  STORE DATA
%  =====================================================================
amp_data = struct();
amp_data.all_VOR            = all_VOR;
amp_data.all_STAT           = all_STAT;
amp_data.temporal_VOR       = all_temporal_VOR;
amp_data.temporal_STAT      = all_temporal_STAT;
amp_data.nasal_VOR          = all_nasal_VOR;
amp_data.nasal_STAT         = all_nasal_STAT;
amp_data.sessions_used      = sessions_used;
amp_data.sessions_skipped   = sessions_skipped;
amp_data.n_sessions         = length(sessions_used);
amp_data.p_ranksum          = p_rs;

%% =====================================================================
%  PLOT + SAVE
%  =====================================================================
plot_amp_figure(amp_data, amp_bin_edges, FontSize, TitleFontSize, TickFontSize, ...
    doShow, doSavePNG, doSaveFIG, doSaveEPS, out_folder);

if doSavePNG || doSaveFIG || doSaveEPS
    save(fullfile(out_folder, 'amplitude_data.mat'), 'amp_data', 'params');
    fprintf('Saved to: %s\n', out_folder);
end

end

%% =====================================================================
%  LOCAL HELPERS
%  =====================================================================

function plot_amp_figure(amp_data, amp_bin_edges, FontSize, TitleFontSize, TickFontSize, ...
        doShow, doSavePNG, doSaveFIG, doSaveEPS, out_folder)

    all_VOR = amp_data.all_VOR;
    all_STAT = amp_data.all_STAT;

    if doShow; fh = figure('Position', [50 50 900 600]);
    else; fh = figure('Position', [50 50 900 600], 'Visible', 'off'); end
    fh.UserData.date_created        = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fh.UserData.generating_function = 'plot_saccade_amplitude_histogram v1.0';

    hold on
    histogram(all_VOR, amp_bin_edges, 'Normalization', 'probability', ...
        'FaceColor', [0.20 0.50 0.80], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    histogram(all_STAT, amp_bin_edges, 'Normalization', 'probability', ...
        'FaceColor', [0.85 0.33 0.10], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    yl = ylim;
    line([median(all_VOR) median(all_VOR)], yl, 'Color', [0.20 0.50 0.80], ...
        'LineWidth', 3, 'LineStyle', '--', 'HandleVisibility', 'off');
    line([median(all_STAT) median(all_STAT)], yl, 'Color', [0.85 0.33 0.10], ...
        'LineWidth', 3, 'LineStyle', '--', 'HandleVisibility', 'off');

    xlabel('Saccade Amplitude (deg)', 'FontSize', FontSize, 'FontWeight', 'bold')
    ylabel('Proportion', 'FontSize', FontSize, 'FontWeight', 'bold')
    set(gca, 'FontSize', TickFontSize)

    legend_vor_str  = sprintf('VOR (n=%d, med=%.1f)', length(all_VOR), median(all_VOR));
    legend_stat_str = sprintf('Spontaneous (n=%d, med=%.1f)', length(all_STAT), median(all_STAT));
    legend({legend_vor_str, legend_stat_str}, 'FontSize', 14, 'Location', 'northeast', 'Box', 'off');

    p_rs = amp_data.p_ranksum;
    if p_rs < 0.001; sig_label = 'p < 0.001';
    elseif p_rs < 0.01; sig_label = 'p < 0.01';
    elseif p_rs < 0.05; sig_label = sprintf('p = %.3f', p_rs);
    else; sig_label = sprintf('p = %.3f n.s.', p_rs);
    end
    text(0.95, 0.72, sig_label, 'Units', 'normalized', 'FontSize', 14, ...
        'HorizontalAlignment', 'right');

    title(sprintf('Saccade Amplitude: VOR vs Spontaneous (%d sessions)', amp_data.n_sessions), ...
        'FontSize', TitleFontSize, 'FontWeight', 'bold')
    box off

    save_fig_trio(fh, fullfile(out_folder, 'saccade_amplitude_histogram'), ...
        doSavePNG, doSaveFIG, doSaveEPS, doShow);
end
function save_fig_trio(fh, filepath_base, doSavePNG, doSaveFIG, doSaveEPS, doShow)
    if doSavePNG; print(fh, [filepath_base '.png'], '-dpng', '-r300'); end
    if doSaveFIG; savefig(fh, [filepath_base '.fig']); end
    if doSaveEPS; print(fh, [filepath_base '.eps'], '-depsc2'); end
    if ~doShow; close(fh); end
end
