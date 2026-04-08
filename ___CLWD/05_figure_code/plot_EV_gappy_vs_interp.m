function [fig_h, params] = plot_EV_gappy_vs_interp(ev_out, varargin)
% plot_EV_gappy_vs_interp  Diagnostic: overlay gappy and interpolated EV traces.
%
%   Shows snippets of the median-filtered slow-phase EV signal (gappy, with
%   physical gaps where saccades/stationary were removed) overlaid with the
%   interpolated uniform-grid signal. Useful for inspecting what the
%   interpolation is doing across gaps.
%
%   INPUTS:
%       ev_out          - output struct from compute_slow_phase_EV
%       varargin        - name-value pairs (see below)
%
%   NAME-VALUE OPTIONS:
%       n_snippets      - number of snippet panels (default: 4)
%       snippet_dur     - duration of each snippet in seconds (default: 5)
%       snippet_starts  - manual start times (overrides auto selection)
%       FontSize        - base font size (default: 18)
%       LineWidth       - line width (default: 2)
%       savePath        - path to save output (default: '', no save)
%       doShow          - display figure (default: 1)
%       doSavePNG       - save PNG (default: 0)
%
%   OUTPUTS:
%       fig_h           - figure handle
%       params          - parameter struct
%
%   USAGE:
%       [ev_out, ~] = compute_slow_phase_EV();
%       plot_EV_gappy_vs_interp(ev_out);
%
%       % With manual snippet start times:
%       plot_EV_gappy_vs_interp(ev_out, 'snippet_starts', [100 200 300 400]);
%
%   DEPENDENCIES:
%       process_varargin
%

%% Defaults
n_snippets      = 4;
snippet_dur     = 5;        % seconds
snippet_starts  = [];       % auto if empty
FontSize        = 18;
LineWidth       = 2;
savePath        = '';
doShow          = 1;
doSavePNG       = 0;

process_varargin(varargin);

%% Store params
params = struct();
params.function_name    = 'plot_EV_gappy_vs_interp';
params.version          = '1.0';
params.date_run         = datestr(now, 'yyyy-mm-dd');
params.n_snippets       = n_snippets;
params.snippet_dur      = snippet_dur;
params.snippet_starts   = snippet_starts;
params.FontSize         = FontSize;

%% Extract data
gappy_tvec  = ev_out.ev_gappy.tvec;
gappy_data  = ev_out.ev_gappy.data;
interp_tvec = ev_out.ev_interp.tvec;
interp_data = ev_out.ev_interp.data;

% Ensure orientations
if size(gappy_tvec, 1) == 1; gappy_tvec = gappy_tvec'; end
if size(gappy_data, 1) ~= 1; gappy_data = gappy_data'; end
if size(interp_tvec, 1) == 1; interp_tvec = interp_tvec'; end
if size(interp_data, 1) ~= 1; interp_data = interp_data'; end

%% Pick snippet start times
t_min = gappy_tvec(1);
t_max = gappy_tvec(end);
total_dur = t_max - t_min;

if isempty(snippet_starts)
    % Space snippets evenly through the recording, avoiding edges
    margin = snippet_dur;
    usable = total_dur - 2 * margin;
    snippet_starts = t_min + margin + linspace(0, usable, n_snippets + 2);
    snippet_starts = snippet_starts(2:end-1);  % trim first and last (edge buffers)
    if length(snippet_starts) > n_snippets
        snippet_starts = snippet_starts(1:n_snippets);
    end
end
n_snippets = length(snippet_starts);

%% Create figure
figW = 1800;
figH = 200 * n_snippets + 100;

if doShow
    fig_h = figure('Position', [50 100 figW figH]);
else
    fig_h = figure('Position', [50 100 figW figH], 'Visible', 'off');
end

for iS = 1:n_snippets
    t_start = snippet_starts(iS);
    t_end   = t_start + snippet_dur;

    ax = subplot(n_snippets, 1, iS); hold on;

    % Gappy trace (discrete samples — plot as individual dots + line)
    g_mask = gappy_tvec >= t_start & gappy_tvec <= t_end;
    g_t = gappy_tvec(g_mask);
    g_d = gappy_data(g_mask);

    % Interpolated trace (continuous uniform grid)
    i_mask = interp_tvec >= t_start & interp_tvec <= t_end;
    i_t = interp_tvec(i_mask);
    i_d = interp_data(i_mask);

    % Plot interpolated first (behind), then gappy on top
    plot(i_t, i_d, '-', 'Color', [0.2 0.6 1.0], 'LineWidth', LineWidth, ...
        'DisplayName', 'Interpolated');
    plot(g_t, g_d, '.', 'Color', [0.8 0.1 0.1], 'MarkerSize', 8, ...
        'DisplayName', 'Gappy (medfilt)');

    xlim([t_start t_end]);
    ylabel('EV (deg/s)', 'FontSize', FontSize);
    set(ax, 'FontSize', FontSize - 2, 'TickDir', 'out', 'Box', 'off');

    if iS == 1
        legend('Interpolated', 'Gappy (medfilt)', 'Location', 'northeast', ...
            'FontSize', FontSize - 4);
        title(sprintf('%s — gappy vs interpolated EV', ev_out.SSN), ...
            'FontSize', FontSize, 'Interpreter', 'none');
    end

    if iS == n_snippets
        xlabel('Time (s)', 'FontSize', FontSize);
    end

    % Mark NaN regions in interpolated signal
    nan_mask_snippet = isnan(i_d);
    if any(nan_mask_snippet)
        % Find contiguous NaN blocks
        nan_diff = diff([0; nan_mask_snippet(:); 0]);
        nan_starts_idx = find(nan_diff == 1);
        nan_ends_idx   = find(nan_diff == -1) - 1;
        for kk = 1:length(nan_starts_idx)
            ns = max(1, nan_starts_idx(kk));
            ne = min(length(i_t), nan_ends_idx(kk));
            patch([i_t(ns) i_t(ne) i_t(ne) i_t(ns)], ...
                  [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)], ...
                  [0.9 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
    end
end

%% Metadata
fig_h.UserData.created_by   = 'plot_EV_gappy_vs_interp';
fig_h.UserData.version      = '1.0';
fig_h.UserData.date_created = datestr(now, 'yyyy-mm-dd');
fig_h.UserData.SSN          = ev_out.SSN;

%% Save
if doSavePNG && ~isempty(savePath)
    timestamp_str = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    save_subfolder = fullfile(savePath, ['EV_gappy_vs_interp_' timestamp_str]);
    if ~exist(save_subfolder, 'dir'); mkdir(save_subfolder); end

    fname = sprintf('EV_gappy_vs_interp_%s.png', ev_out.SSN);
    print(fig_h, fullfile(save_subfolder, fname), '-dpng', '-r300');
    fprintf('Saved: %s\n', fullfile(save_subfolder, fname));

    save(fullfile(save_subfolder, 'EV_gappy_vs_interp_params.mat'), 'params');
end

if ~doShow
    close(fig_h);
end

end

%% Version History
% 2026-03-18  v1.0 Initial version.
