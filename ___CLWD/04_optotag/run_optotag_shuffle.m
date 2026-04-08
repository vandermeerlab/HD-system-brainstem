%% run_optotag_shuffle.m
% Runs compute_optotag_shuffle on all 7 optotagged cells and saves results.
% Requires cell_data in workspace (from plot_optotag_laser_PETHs v2.3+).

%% Config
cfg_shuff = [];
cfg_shuff.nShuffles = 1000;
cfg_shuff.verbose = 0;

nCells = length(cell_data);
shuff_results = cell(nCells, 1);

%% Run shuffle on each cell
for iC = 1:nCells
    fprintf('Cell %d/%d...\n', iC, nCells);
    [shuff_results{iC}, ~] = compute_optotag_shuffle( ...
        cell_data{iC}.spike_times, ...
        cell_data{iC}.laser_on, ...
        cell_data{iC}.session_duration, ...
        cfg_shuff);
    fprintf('  percentile = %.4f\n', shuff_results{iC}.percentile_overall);
end

%% Print summary table
fprintf('\nCell\tPeak FR (real)\tPeak FR (shuff mean)\tPercentile\n');
fprintf('----\t--------------\t--------------------\t----------\n');
for iC = 1:nCells
    fprintf('%d\t%.2f\t\t%.2f\t\t\t%.4f\n', iC, ...
        shuff_results{iC}.peak_FR_real, ...
        mean(shuff_results{iC}.peak_FR_shuff), ...
        shuff_results{iC}.percentile_overall);
end

%% Save
save_dir = 'C:\Jeff\U01\datatouse\___analysis_output\figS4_optotag_PETHs_2026-03-18';
if ~exist(save_dir, 'dir'); mkdir(save_dir); end
save_path = fullfile(save_dir, 'figS4_shuffle_results.mat');
save(save_path, 'shuff_results', 'cfg_shuff', '-v7.3');
fprintf('\nSaved to %s\n', save_path);
