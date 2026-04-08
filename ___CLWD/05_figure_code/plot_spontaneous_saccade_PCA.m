function [pca_results, params] = plot_spontaneous_saccade_PCA(varargin)
% plot_spontaneous_saccade_PCA  v1.0
% PCA of spontaneous saccade PETHs showing bidirectional modulation.
% Analog of Figure S5B but restricted to platform-stationary saccades.
%
%   Takes preferred and non-preferred spontaneous saccade PETHs,
%   runs PCA on the stacked matrix, and plots PC1 coefficients
%   for preferred vs non-preferred. Negative slope demonstrates
%   opposite modulation across the population.
%
%   INPUTS (name-value pairs via varargin):
%     Z                 - Z struct from plot_stationary_saccade_heatmap
%                         (if empty, loads from population_PETH_data.mat)
%     pop_data_path     - path to population_PETH_data.mat (fallback)
%     n_example_insets  - how many example PETHs to show as insets (default: 4)
%     doSavePNG/FIG/EPS - save flags (default: 1)
%     doShow            - display figure (default: 1)
%     save_dir          - output root
%
%   OUTPUTS:
%     pca_results - struct with PC coefficients, explained variance, fit stats
%     params      - struct of parameter values
%
%   DEPENDENCIES:
%     process_varargin

%% =====================================================================
%  DEFAULTS
%  =====================================================================
Z                 = [];
pop_data_path     = '';
n_example_insets  = 0;
save_dir          = 'C:\Jeff\U01\datatouse\___analysis_output';

% Figure
FontSize          = 24;
TitleFontSize     = 22;
TickFontSize      = 16;
MarkerSize_scatter = 60;
color_pref        = [0.85 0.15 0.15];   % red
color_nonpref     = [0.15 0.65 0.15];   % green
color_fit_line    = [0.3 0.3 0.3];

% Output control
doSavePNG         = 1;
doSaveFIG         = 1;
doSaveEPS         = 1;
doShow            = 1;

process_varargin(varargin);

% Manual extraction for struct/string params that process_varargin can miss
for iV = 1:2:length(varargin)
    if ischar(varargin{iV})
        switch varargin{iV}
            case 'Z';              Z = varargin{iV+1};
            case 'pop_data_path';  pop_data_path = varargin{iV+1};
            case 'precomputed_Z';  Z = varargin{iV+1};  % alias
        end
    end
end

%% =====================================================================
%  PARAMS
%  =====================================================================
params = struct();
params.function_name    = 'plot_spontaneous_saccade_PCA';
params.version          = '1.3';
params.date_run         = datestr(now, 'yyyy-mm-dd HH:MM:SS');
params.n_example_insets = n_example_insets;

%% =====================================================================
%  GET DATA
%  =====================================================================
has_Z = ~isempty(Z) && isstruct(Z) && isfield(Z, 'FR_preferred_smooth');
has_pop = ~isempty(pop_data_path) && ischar(pop_data_path) && length(pop_data_path) > 0;

if has_Z
    % Use Z struct from heatmap function
    pref_mat    = Z.FR_preferred_smooth;
    nonpref_mat = Z.FR_nonpreferred_smooth;
    peth_centers = Z.binCenters;
    n_neurons    = Z.n_neurons;
    n_sessions   = Z.n_sessions;
    min_saccades = Z.min_saccades;
    neuron_names = Z.neuron_names;
    data_source  = 'Z struct from plot_stationary_saccade_heatmap';
elseif has_pop
    if ~exist(pop_data_path, 'file')
        error('pop_data_path file not found: %s', pop_data_path);
    end
    % Load from population_PETH_data.mat
    loaded = load(pop_data_path, 'pop_data');
    pd = loaded.pop_data;
    
    % Ensure neurons x bins orientation
    % peth_centers should be 1 x nBins
    peth_centers = pd.peth_centers(:)';
    nBins_expected = length(peth_centers);
    
    % pref_spont_raw: if size is (nBins x nNeurons), transpose to (nNeurons x nBins)
    raw_pref = pd.pref_spont_raw;
    if size(raw_pref, 1) == nBins_expected && size(raw_pref, 2) ~= nBins_expected
        pref_mat = raw_pref';
    elseif size(raw_pref, 2) == nBins_expected
        pref_mat = raw_pref;
    else
        error('Unexpected pref_spont_raw dimensions: [%d x %d], expected %d bins', ...
            size(raw_pref,1), size(raw_pref,2), nBins_expected);
    end
    
    raw_nonpref = pd.nonpref_spont_raw;
    if size(raw_nonpref, 1) == nBins_expected && size(raw_nonpref, 2) ~= nBins_expected
        nonpref_mat = raw_nonpref';
    elseif size(raw_nonpref, 2) == nBins_expected
        nonpref_mat = raw_nonpref;
    else
        error('Unexpected nonpref_spont_raw dimensions: [%d x %d], expected %d bins', ...
            size(raw_nonpref,1), size(raw_nonpref,2), nBins_expected);
    end
    peth_centers = pd.peth_centers(:)';
    n_neurons    = pd.n_neurons;
    n_sessions   = pd.n_sessions;
    min_saccades = 10;
    % Extract names
    neuron_names = cell(n_neurons, 1);
    for ii = 1:n_neurons
        neuron_names{ii} = pd.names{ii};
    end
    data_source = pop_data_path;
else
    error('Provide either Z struct or pop_data_path. pop_data_path = ''%s''', pop_data_path);
end

params.data_source = data_source;
params.n_neurons   = n_neurons;

fprintf('PCA on %d neurons (%d sessions)\n', n_neurons, n_sessions);

%% =====================================================================
%  OUTPUT FOLDER
%  =====================================================================
run_date_str = datestr(now, 'yyyy-mm-dd');
run_num = 1;
out_folder = fullfile(save_dir, sprintf('spont_PCA_%s_run%d', run_date_str, run_num));
while exist(out_folder, 'dir')
    run_num = run_num + 1;
    out_folder = fullfile(save_dir, sprintf('spont_PCA_%s_run%d', run_date_str, run_num));
end
if doSavePNG || doSaveFIG || doSaveEPS; mkdir(out_folder); end
params.output_folder = out_folder;

%% =====================================================================
%  RUN PCA — matching Matt's approach from SESSION_ahv_glmfit
%  =====================================================================
% Step 1: Normalize each neuron's PETH by its max FR
%          (like NnormSmooth/TnormSmooth from makeSaccadeHeatPlot)
pref_max  = max(pref_mat, [], 2);
pref_max(pref_max == 0) = 1;  % avoid divide by zero
pref_norm = pref_mat ./ pref_max;

nonpref_max  = max(nonpref_mat, [], 2);
nonpref_max(nonpref_max == 0) = 1;
nonpref_norm = nonpref_mat ./ nonpref_max;

% Step 2: Smooth (matching makeSaccadeHeatPlot's smoothdata call)
pref_norm_smooth    = smoothdata(pref_norm, 2);
nonpref_norm_smooth = smoothdata(nonpref_norm, 2);

% Step 3: Stack as nBins x 2*nNeurons (Matt's convention)
%   x1 = NnormSmooth'; x2 = TnormSmooth'; x = cat(2, x1, x2);
%   rows = time bins (observations), columns = neurons (variables)
x_stacked = cat(2, pref_norm_smooth', nonpref_norm_smooth');

% Remove any columns with NaN
valid_cols = ~any(isnan(x_stacked), 1);
x_clean = x_stacked(:, valid_cols);

% Track which neurons survived NaN removal
n_total_cols = 2 * n_neurons;
valid_pref_idx    = valid_cols(1:n_neurons);
valid_nonpref_idx = valid_cols(n_neurons+1:n_total_cols);

% Step 4: PCA — coeff columns are the loadings per neuron
[coeff_pca, score_pca, latent_pca, ~, explained_pca] = pca(x_clean);

% Orient PC1 so that the saccade burst (peak near t=0) is positive
[~, t0_bin] = min(abs(peth_centers));  % bin closest to saccade onset
if score_pca(t0_bin, 1) < 0
    coeff_pca(:, 1) = -coeff_pca(:, 1);
    score_pca(:, 1) = -score_pca(:, 1);
end

% Step 5: Read loadings directly from coeff(:,1)
%   coeff is (2*nValid) x nComponents
%   First half = preferred side loadings, second half = non-preferred
n_valid_pref    = sum(valid_pref_idx);
n_valid_nonpref = sum(valid_nonpref_idx);

pc1_pref    = nan(n_neurons, 1);
pc1_nonpref = nan(n_neurons, 1);
pc1_pref(valid_pref_idx)       = coeff_pca(1:n_valid_pref, 1);
pc1_nonpref(valid_nonpref_idx) = coeff_pca(n_valid_pref+1:n_valid_pref+n_valid_nonpref, 1);

% Linear fit
valid_both = ~isnan(pc1_pref) & ~isnan(pc1_nonpref);
fit_coeffs = polyfit(pc1_pref(valid_both), pc1_nonpref(valid_both), 1);
fit_r = corr(pc1_pref(valid_both), pc1_nonpref(valid_both));

fprintf('  PC1 explains %.1f%% of variance\n', explained_pca(1));
fprintf('  Slope of fit: %.3f\n', fit_coeffs(1));
fprintf('  Correlation:  r = %.3f\n', fit_r);

%% =====================================================================
%  STORE RESULTS
%  =====================================================================
pca_results = struct();
pca_results.pc1_pref         = pc1_pref;
pca_results.pc1_nonpref      = pc1_nonpref;
pca_results.pc1_timecourse   = score_pca(:, 1);
pca_results.coeff            = coeff_pca;
pca_results.latent           = latent_pca;
pca_results.explained        = explained_pca;
pca_results.fit_slope        = fit_coeffs(1);
pca_results.fit_intercept    = fit_coeffs(2);
pca_results.fit_r            = fit_r;
pca_results.peth_centers     = peth_centers;
pca_results.neuron_names     = neuron_names;
pca_results.n_neurons        = n_neurons;
pca_results.pref_norm_smooth = pref_norm_smooth;
pca_results.nonpref_norm_smooth = nonpref_norm_smooth;

%% =====================================================================
%  PLOT
%  =====================================================================
if doShow
    fh = figure('Position', [50 50 900 750]);
else
    fh = figure('Position', [50 50 900 750], 'Visible', 'off');
end
fh.UserData.date_created        = datestr(now, 'yyyy-mm-dd HH:MM:SS');
fh.UserData.generating_function = 'plot_spontaneous_saccade_PCA v1.3';

hold on

% Scatter
scatter(pc1_pref, pc1_nonpref, MarkerSize_scatter, [0.4 0.4 0.4], 'filled', ...
    'MarkerFaceAlpha', 0.6);

% Fit line
xl = [min(pc1_pref)*1.1 max(pc1_pref)*1.1];
fit_y = polyval(fit_coeffs, xl);
plot(xl, fit_y, '-', 'Color', color_fit_line, 'LineWidth', 3);

% Zero lines
xline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'HandleVisibility', 'off');
yline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Preferred saccade PC1 coeff', 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', color_pref)
ylabel('Non-preferred saccade PC1 coeff', 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', color_nonpref)
set(gca, 'FontSize', TickFontSize)
box off

% Annotations
text(0.05, 0.95, sprintf('r = %.3f', fit_r), 'Units', 'normalized', ...
    'FontSize', 20, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
text(0.05, 0.88, sprintf('slope = %.2f', fit_coeffs(1)), 'Units', 'normalized', ...
    'FontSize', 18, 'VerticalAlignment', 'top');
text(0.05, 0.81, sprintf('PC1: %.0f%% var', explained_pca(1)), 'Units', 'normalized', ...
    'FontSize', 18, 'VerticalAlignment', 'top');

title(sprintf('Spontaneous Saccade PCA (%d neurons, >= %d saccades/dir)', ...
    n_neurons, min_saccades), 'FontSize', TitleFontSize, 'FontWeight', 'bold');

% NOTE: Example insets disabled. Re-enable with n_example_insets > 0.
% Currently set to 0 to avoid cluttered display.

%% =====================================================================
%  SAVE
%  =====================================================================
if doSavePNG || doSaveFIG || doSaveEPS
    save_fig_trio(fh, fullfile(out_folder, 'spontaneous_saccade_PCA'), ...
        doSavePNG, doSaveFIG, doSaveEPS, doShow);
    save(fullfile(out_folder, 'pca_results.mat'), 'pca_results', 'params');
    fprintf('Saved to: %s\n', out_folder);
end

end

%% =====================================================================
%  LOCAL HELPERS
%  =====================================================================
function save_fig_trio(fh, fp, doP, doF, doE, doS)
    if doP; print(fh, [fp '.png'], '-dpng', '-r300'); end
    if doF; savefig(fh, [fp '.fig']); end
    if doE; print(fh, [fp '.eps'], '-depsc2'); end
    if ~doS; close(fh); end
end
