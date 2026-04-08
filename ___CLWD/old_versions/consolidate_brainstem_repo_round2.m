%% Brainstem Repo Consolidation — Round 2
%  Generated: 2026-04-06 (Session 7)
%
%  Creates umbrella folders for AHV, eye_movements, stimuli.
%  Moves domain-specific plotting into domain folders.
%  Slims plotting/ to general-purpose only.
%
%  SET dryRun = false TO EXECUTE. Default is true (print-only).

dryRun = false;  % <<< FLIP TO false WHEN READY

repo = 'C:\Jeff\Github_repositories\HD-system-brainstem';

% --- New umbrella folders and subfolders ---
newFolders = {
    fullfile(repo, 'AHV')
    fullfile(repo, 'AHV', 'graham_classification')
    fullfile(repo, 'AHV', 'regression')
    fullfile(repo, 'AHV', 'tuning_curves')
    fullfile(repo, 'AHV', 'glm')
    fullfile(repo, 'eye_movements')
    fullfile(repo, 'eye_movements', 'slow_phase')
    fullfile(repo, 'eye_movements', 'saccades')
    fullfile(repo, 'eye_movements', 'pupil_processing')
    fullfile(repo, 'stimuli')
};

for iF = 1:length(newFolders)
    if ~exist(newFolders{iF}, 'dir')
        if dryRun
            fprintf('[DRY RUN] MKDIR: %s\n', newFolders{iF});
        else
            mkdir(newFolders{iF});
            fprintf('CREATED: %s\n', newFolders{iF});
        end
    end
end

% =====================================================================
%  MOVE LIST
%  Each entry: {source_relative_to_repo, destination_folder}
%  Source can be a folder (merges all .m files) or a single file.
% =====================================================================

moves = {};

% -----------------------------------------------------------------
%  1. AHV UMBRELLA — merge 4 existing folders into subfolders
% -----------------------------------------------------------------

% AHV graham defs/ → AHV/graham_classification/
moves{end+1,1} = 'AHV graham defs';
moves{end,  2} = fullfile(repo, 'AHV', 'graham_classification');

% AHV regression/ → AHV/regression/
moves{end+1,1} = 'AHV regression';
moves{end,  2} = fullfile(repo, 'AHV', 'regression');

% AHV tuning curves/ → AHV/tuning_curves/
moves{end+1,1} = 'AHV tuning curves';
moves{end,  2} = fullfile(repo, 'AHV', 'tuning_curves');

% ahv_glmfit/ → AHV/glm/
moves{end+1,1} = 'ahv_glmfit';
moves{end,  2} = fullfile(repo, 'AHV', 'glm');

% AHV-specific plotting → AHV/ (loose at umbrella level)
moves{end+1,1} = fullfile('plotting', 'plotAHVscatter.m');
moves{end,  2} = fullfile(repo, 'AHV');
moves{end+1,1} = fullfile('plotting', 'plot_AHVtuning.m');
moves{end,  2} = fullfile(repo, 'AHV', 'tuning_curves');
moves{end+1,1} = fullfile('plotting', 'plotAHVtuningplussaccadePETHs.m');
moves{end,  2} = fullfile(repo, 'AHV');
moves{end+1,1} = fullfile('plotting', 'plotTCbyMaxBin.m');
moves{end,  2} = fullfile(repo, 'AHV');

% -----------------------------------------------------------------
%  2. EYE_MOVEMENTS UMBRELLA — merge 4 existing folders + plotting
% -----------------------------------------------------------------

% slow phase/ → eye_movements/slow_phase/
moves{end+1,1} = 'slow phase';
moves{end,  2} = fullfile(repo, 'eye_movements', 'slow_phase');

% saccade peths/ → eye_movements/saccades/
moves{end+1,1} = 'saccade peths';
moves{end,  2} = fullfile(repo, 'eye_movements', 'saccades');

% saccade quantification/ → eye_movements/saccades/
moves{end+1,1} = 'saccade quantification';
moves{end,  2} = fullfile(repo, 'eye_movements', 'saccades');

% pupil processing/ → eye_movements/pupil_processing/
moves{end+1,1} = 'pupil processing';
moves{end,  2} = fullfile(repo, 'eye_movements', 'pupil_processing');

% Saccade/eye-specific plotting → eye_movements/saccades/
moves{end+1,1} = fullfile('plotting', 'plotSaccadeTrace.m');
moves{end,  2} = fullfile(repo, 'eye_movements', 'saccades');
moves{end+1,1} = fullfile('plotting', 'plotSaccades.m');
moves{end,  2} = fullfile(repo, 'eye_movements', 'saccades');
moves{end+1,1} = fullfile('plotting', 'plotSaccadesAndAHVsinglePlot5.m');
moves{end,  2} = fullfile(repo, 'eye_movements', 'saccades');
moves{end+1,1} = fullfile('plotting', 'SpikePETHsaccades.m');
moves{end,  2} = fullfile(repo, 'eye_movements', 'saccades');
moves{end+1,1} = fullfile('plotting', 'saccadePETHsubplots.m');
moves{end,  2} = fullfile(repo, 'eye_movements', 'saccades');

% EV plotting → eye_movements/slow_phase/
moves{end+1,1} = fullfile('plotting', 'plot_EV_TC_from_workspace.m');
moves{end,  2} = fullfile(repo, 'eye_movements', 'slow_phase');

% DTN analyses/ (1 file) → eye_movements/saccades/
moves{end+1,1} = 'DTN analyses';
moves{end,  2} = fullfile(repo, 'eye_movements', 'saccades');

% -----------------------------------------------------------------
%  3. STIMULI UMBRELLA — merge 3 small folders
% -----------------------------------------------------------------

% airpuff/ → stimuli/
moves{end+1,1} = 'airpuff';
moves{end,  2} = fullfile(repo, 'stimuli');

% audio and laser stimuli/ → stimuli/
moves{end+1,1} = 'audio and laser stimuli';
moves{end,  2} = fullfile(repo, 'stimuli');

% OKR/ → stimuli/
moves{end+1,1} = 'OKR';
moves{end,  2} = fullfile(repo, 'stimuli');

% -----------------------------------------------------------------
%  4. SLIM PLOTTING/ — move non-general files out
% -----------------------------------------------------------------

% brainstem_megaplot_allsessions.m → batch_runners/
moves{end+1,1} = fullfile('plotting', 'brainstem_megaplot_allsessions.m');
moves{end,  2} = fullfile(repo, 'batch_runners');

% subplot_solution.m → for_review/
moves{end+1,1} = fullfile('plotting', 'subplot_solution.m');
moves{end,  2} = fullfile(repo, 'for_review');

% -----------------------------------------------------------------
%  5. SMALL MERGES
% -----------------------------------------------------------------

% cell_selection/ (3 files) → checks and summary numbers/
moves{end+1,1} = 'cell_selection';
moves{end,  2} = fullfile(repo, 'checks and summary numbers');


% =====================================================================
%  EXECUTE (or dry-run)
% =====================================================================

fprintf('\n========================================\n');
if dryRun
    fprintf('  DRY RUN — no files will be moved\n');
else
    fprintf('  EXECUTING — files will be moved!\n');
end
fprintf('========================================\n\n');

nMoved = 0;
nSkipped = 0;
nFailed = 0;
folderMerges = {};

for iM = 1:size(moves, 1)
    srcRel = moves{iM, 1};
    dstDir = moves{iM, 2};
    srcFull = fullfile(repo, srcRel);

    % --- Handle folder merges (move all .m files from folder) ---
    if exist(srcFull, 'dir') == 7
        % Grab ALL files (not just .m) to preserve .mat, .txt, etc.
        folderFiles = dir(srcFull);
        folderFiles = folderFiles(~ismember({folderFiles.name}, {'.', '..'}));
        % Separate files from subdirectories
        isSubdir = [folderFiles.isdir];
        fileEntries = folderFiles(~isSubdir);

        if isempty(fileEntries)
            fprintf('  EMPTY FOLDER (no files): %s\n', srcRel);
        else
            fprintf('\n--- Merging folder: %s → %s ---\n', srcRel, dstDir);
            for jj = 1:length(fileEntries)
                srcFile = fullfile(srcFull, fileEntries(jj).name);
                dstFile = fullfile(dstDir, fileEntries(jj).name);
                if exist(dstFile, 'file')
                    fprintf('  CONFLICT (already exists): %s\n', fileEntries(jj).name);
                    nSkipped = nSkipped + 1;
                    continue;
                end
                if dryRun
                    fprintf('  [DRY] %s\n', fileEntries(jj).name);
                else
                    try
                        movefile(srcFile, dstFile);
                        fprintf('  MOVED: %s\n', fileEntries(jj).name);
                    catch ME
                        fprintf('  FAILED: %s — %s\n', fileEntries(jj).name, ME.message);
                        nFailed = nFailed + 1;
                    end
                end
                nMoved = nMoved + 1;
            end
        end

        % Check for subdirectories to move as whole folders
        subdirEntries = folderFiles(isSubdir);
        for jj = 1:length(subdirEntries)
            srcSubdir = fullfile(srcFull, subdirEntries(jj).name);
            dstSubdir = fullfile(dstDir, subdirEntries(jj).name);
            if exist(dstSubdir, 'dir')
                fprintf('  CONFLICT (subdir already exists): %s\n', subdirEntries(jj).name);
                nSkipped = nSkipped + 1;
                continue;
            end
            if dryRun
                fprintf('  [DRY] SUBDIR: %s/\n', subdirEntries(jj).name);
            else
                try
                    movefile(srcSubdir, dstSubdir);
                    fprintf('  MOVED SUBDIR: %s/\n', subdirEntries(jj).name);
                catch ME
                    fprintf('  FAILED SUBDIR: %s — %s\n', subdirEntries(jj).name, ME.message);
                    nFailed = nFailed + 1;
                end
            end
            nMoved = nMoved + 1;
        end

        folderMerges{end+1} = srcFull; %#ok<SAGROW>
        continue;
    end

    % --- Handle single file move ---
    if ~exist(srcFull, 'file')
        fprintf('  SKIP (not found): %s\n', srcRel);
        nSkipped = nSkipped + 1;
        continue;
    end

    [~, fname, fext] = fileparts(srcFull);
    dstFile = fullfile(dstDir, [fname fext]);

    if exist(dstFile, 'file')
        fprintf('  CONFLICT (already exists at dest): %s\n', srcRel);
        nSkipped = nSkipped + 1;
        continue;
    end

    if dryRun
        fprintf('[DRY] MOVE: %s\n       → %s\n', srcRel, dstDir);
    else
        try
            movefile(srcFull, dstFile);
            fprintf('MOVED: %s → %s\n', srcRel, dstDir);
        catch ME
            fprintf('FAILED: %s — %s\n', srcRel, ME.message);
            nFailed = nFailed + 1;
        end
    end
    nMoved = nMoved + 1;
end

% --- Summary ---
fprintf('\n========================================\n');
fprintf('  SUMMARY\n');
fprintf('========================================\n');
if dryRun
    fprintf('  Would move: %d items\n', nMoved);
else
    fprintf('  Moved: %d items\n', nMoved);
end
fprintf('  Skipped (not found or conflict): %d\n', nSkipped);
if ~dryRun
    fprintf('  Failed: %d\n', nFailed);
end

% --- Report empty source folders ---
if ~isempty(folderMerges)
    fprintf('\n  Source folders to delete if empty:\n');
    for ii = 1:length(folderMerges)
        remaining = dir(folderMerges{ii});
        remaining = remaining(~ismember({remaining.name}, {'.', '..'}));
        if isempty(remaining)
            fprintf('    [EMPTY — safe to delete] %s\n', folderMerges{ii});
        else
            fprintf('    [%d items remain] %s\n', length(remaining), folderMerges{ii});
        end
    end
end

% --- What stays in plotting/ ---
fprintf('\n  Files remaining in plotting/ after consolidation:\n');
pltDir = fullfile(repo, 'plotting');
pltRemaining = dir(fullfile(pltDir, '*.m'));
movedNames = {};
for iM = 1:size(moves,1)
    [~, fn, fe] = fileparts(moves{iM,1});
    movedNames{end+1} = [fn fe]; %#ok<SAGROW>
end
for jj = 1:length(pltRemaining)
    if ~ismember(pltRemaining(jj).name, movedNames)
        fprintf('    %s\n', pltRemaining(jj).name);
    end
end

% --- New folder structure overview ---
fprintf('\n  New umbrella structure:\n');
umbrellas = {'AHV', 'eye_movements', 'stimuli'};
for ii = 1:length(umbrellas)
    uDir = fullfile(repo, umbrellas{ii});
    fprintf('    %s/\n', umbrellas{ii});
    if exist(uDir, 'dir')
        subs = dir(uDir);
        subs = subs([subs.isdir] & ~ismember({subs.name}, {'.', '..'}));
        for jj = 1:length(subs)
            fprintf('      %s/\n', subs(jj).name);
        end
    end
end

fprintf('\n');
if dryRun
    fprintf('>>> Set dryRun = false and re-run to execute. <<<\n');
else
    fprintf('>>> Done! Delete empty source folders, then run CLWD_inventory_ver_1_0 to verify. <<<\n');
end
