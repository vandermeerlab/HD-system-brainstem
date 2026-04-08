%% Brainstem Repo Reorganization Script
%  Generated: 2026-04-06 (Session 7)
%  
%  This script reorganizes HD-system-brainstem by:
%    1. Archiving old numbered versions → old_versions_archive/
%    2. Moving ambiguous files → for_review/
%    3. Splitting helper functions/ into focused subfolders
%    4. Merging small folders into related larger ones
%    5. Cleaning up root-level loose files
%
%  SET dryRun = false TO EXECUTE. Default is true (print-only).

dryRun = false;  % <<< EXECUTING FOR REAL

repo = 'C:\Jeff\Github_repositories\HD-system-brainstem';
archiveDir   = fullfile(repo, 'old_versions_archive');
reviewDir    = fullfile(repo, 'for_review');
sessionUtils = fullfile(repo, 'session_utils');
batchRunners = fullfile(repo, 'batch_runners');
dataFixes    = fullfile(repo, 'data_fixes');

% --- Create destination folders ---
newFolders = {archiveDir, reviewDir, sessionUtils, batchRunners, dataFixes};
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
%  MOVE LIST: {source, destination}
%  Source paths are relative to repo root.
%  Destination is a folder — file keeps its name.
% =====================================================================

moves = {};

% -----------------------------------------------------------------
%  1. ARCHIVE — old numbered versions
% -----------------------------------------------------------------

% plotting/ megaplot iterations (brainstem_megaplot.m is the current one)
for ii = 1:9
    moves{end+1,1} = fullfile('plotting', sprintf('plot_AHVtuningMegaplot%d.m', ii));
    moves{end,  2} = archiveDir;
end
% The unnumbered original
moves{end+1,1} = fullfile('plotting', 'plot_AHVtuningMegaplot.m');
moves{end,  2} = archiveDir;

% megaplot10 is equivalent to brainstem_megaplot — archive it too
moves{end+1,1} = fullfile('plotting', 'plot_AHVtuningMegaplot10.m');
moves{end,  2} = archiveDir;

% plotSaccadesAndAHVsinglePlot iterations (keep 5)
for ii = 1:4
    if ii == 1
        moves{end+1,1} = fullfile('plotting', 'plotSaccadesAndAHVsinglePlot.m');
    else
        moves{end+1,1} = fullfile('plotting', sprintf('plotSaccadesAndAHVsinglePlot%d.m', ii));
    end
    moves{end,2} = archiveDir;
end

% SpikePETH older variants (SpikePETH_either.m is current)
moves{end+1,1} = fullfile('plotting', 'SpikePETH_either2.m');
moves{end,  2} = archiveDir;
moves{end+1,1} = fullfile('plotting', 'SpikePETHvdm.m');
moves{end,  2} = archiveDir;

% pupil processing/ — detectSaccadesManualCheck 2-6 (keep 7)
for ii = 2:6
    moves{end+1,1} = fullfile('pupil processing', sprintf('detectSaccadesManualCheck%d.m', ii));
    moves{end,  2} = archiveDir;
end
moves{end+1,1} = fullfile('pupil processing', 'calibrate_eye_position_single_session_ver1.m');
moves{end,  2} = archiveDir;

% slow phase/ — old versions (keep getSlowPhaseData2, getSlowPhaseTC3)
moves{end+1,1} = fullfile('slow phase', 'getSlowPhaseData.m');
moves{end,  2} = archiveDir;
moves{end+1,1} = fullfile('slow phase', 'getSlowPhaseData_old.m');
moves{end,  2} = archiveDir;
moves{end+1,1} = fullfile('slow phase', 'getSlowPhaseTC.m');
moves{end,  2} = archiveDir;
moves{end+1,1} = fullfile('slow phase', 'getSlowPhaseTC2.m');
moves{end,  2} = archiveDir;

% helper functions/ — superseded loaders & event label versions
moves{end+1,1} = fullfile('helper functions', 'SortBrainstemEventLabels.m');
moves{end,  2} = archiveDir;
moves{end+1,1} = fullfile('helper functions', 'SortBrainstemEventLabels2.m');
moves{end,  2} = archiveDir;
moves{end+1,1} = fullfile('helper functions', 'SortBrainstemEventLabelsOld.m');
moves{end,  2} = archiveDir;
moves{end+1,1} = fullfile('helper functions', 'LoadSpikesJeff.m');
moves{end,  2} = archiveDir;
moves{end+1,1} = fullfile('helper functions', 'LoadSpikesJeff2.m');
moves{end,  2} = archiveDir;
moves{end+1,1} = fullfile('helper functions', 'HF_CreateBlankRecordingKeys.m');
moves{end,  2} = archiveDir;
moves{end+1,1} = fullfile('helper functions', 'find_tt_nums.m');
moves{end,  2} = archiveDir;

% other versions of shadederrorbar/ — entire folder contents
moves{end+1,1} = fullfile('other versions of shadederrorbar');
moves{end,  2} = archiveDir;  % handled specially below as folder

% old code versions/ — entire folder (12 superseded saccade sig tests)
moves{end+1,1} = fullfile('old code versions');
moves{end,  2} = archiveDir;  % handled specially below as folder

% -----------------------------------------------------------------
%  2. FOR_REVIEW — ambiguous files, not typically used
% -----------------------------------------------------------------
moves{end+1,1} = fullfile('helper functions', 'MultiRasterJeff.m');
moves{end,  2} = reviewDir;
moves{end+1,1} = fullfile('helper functions', 'PlotSpikeRaster2_Jeff.m');
moves{end,  2} = reviewDir;
moves{end+1,1} = fullfile('helper functions', 'CheckHFkeyss.m');
moves{end,  2} = reviewDir;
moves{end+1,1} = fullfile('helper functions', 'simplePETH.m');
moves{end,  2} = reviewDir;
moves{end+1,1} = fullfile('helper functions', 'findSessionStartTime.m');
moves{end,  2} = reviewDir;
moves{end+1,1} = fullfile('helper functions', 'TuningCurvesSDF.m');
moves{end,  2} = reviewDir;

% -----------------------------------------------------------------
%  3. SPLIT helper functions/ into focused subfolders
% -----------------------------------------------------------------

% --- session_utils/ ---
sessionUtilFiles = {
    'process_varargin.m'
    'HD_GetSSN.m'
    'EvalKeys.m'
    'getDataPath.m'
    'sortSessionsByDate.m'
    'findSessionEndTime.m'
    'CountNumNeuronsPerSession.m'
    'find_tt_nums2.m'
    'tally_hemisphere_tfilelist.m'
    'isOdd.m'
    'get_keys_list.m'
};
for ii = 1:length(sessionUtilFiles)
    moves{end+1,1} = fullfile('helper functions', sessionUtilFiles{ii});
    moves{end,  2} = sessionUtils;
end

% --- batch_runners/ ---
batchFiles = {
    'countSaccadesAllSessions.m'
    'clusterQual_allsessions.m'
    'collect_all_cluQual.m'
    'getEventLabels_Allsessions.m'
    'getLaserTimesAllSessions.m'
    'run_laserEvents_all_sessions.m'
    'findSteadyPeriodsManual_AllSessions.m'
    'sessionLoopTemplate.m'
};
for ii = 1:length(batchFiles)
    moves{end+1,1} = fullfile('helper functions', batchFiles{ii});
    moves{end,  2} = batchRunners;
end

% --- data_fixes/ ---
dataFixFiles = {
    'correctSaccadeNaNs.m'
    'multiply_by_50.m'
    'HF_changefilenames.m'
    'deletemp4files.m'
};
for ii = 1:length(dataFixFiles)
    moves{end+1,1} = fullfile('helper functions', dataFixFiles{ii});
    moves{end,  2} = dataFixes;
end

% --- Move domain-specific helpers to their natural homes ---
% Plotting-related
moves{end+1,1} = fullfile('helper functions', 'plotSTART_STOP_lines.m');
moves{end,  2} = fullfile(repo, 'plotting');
moves{end+1,1} = fullfile('helper functions', 'generate_EYE_AHV_raster.m');
moves{end,  2} = fullfile(repo, 'plotting');
moves{end+1,1} = fullfile('helper functions', 'moveLabel.m');
moves{end,  2} = fullfile(repo, 'plotting');
moves{end+1,1} = fullfile('helper functions', 'ylabel_spacing.m');
moves{end,  2} = fullfile(repo, 'plotting');

% Saccade-related
moves{end+1,1} = fullfile('helper functions', 'getSaccadeAmplitudeInfo.m');
moves{end,  2} = fullfile(repo, 'saccade quantification');

% Eye tracking diagnostic
moves{end+1,1} = fullfile('helper functions', 'show_ipsi_contra_nasal_temporal_examples.m');
moves{end,  2} = fullfile(repo, 'pupil processing');

% -----------------------------------------------------------------
%  4. MERGE small folders into related larger ones
% -----------------------------------------------------------------

% AHV score/ (2 files) → AHV graham defs/
moves{end+1,1} = fullfile('AHV score');
moves{end,  2} = fullfile(repo, 'AHV graham defs');

% cell summary statistics/ (1 file) → checks and summary numbers/
moves{end+1,1} = fullfile('cell summary statistics');
moves{end,  2} = fullfile(repo, 'checks and summary numbers');

% util/ (2 files) → helper functions/
moves{end+1,1} = fullfile('util');
moves{end,  2} = fullfile(repo, 'helper functions');

% postprocessing functions/ (1 file) → pupil processing/
moves{end+1,1} = fullfile('postprocessing functions');
moves{end,  2} = fullfile(repo, 'pupil processing');

% -----------------------------------------------------------------
%  5. ROOT-LEVEL cleanup
% -----------------------------------------------------------------

% temp.m — 0 bytes, delete
moves{end+1,1} = 'temp.m';
moves{end,  2} = '__DELETE__';

% fd_to_fdr.m at root — duplicate of helper functions/fdtofdr.m
moves{end+1,1} = 'fd_to_fdr.m';
moves{end,  2} = archiveDir;

% saccade_peth_MVDM.m — belongs with saccade code
moves{end+1,1} = 'saccade_peth_MVDM.m';
moves{end,  2} = fullfile(repo, 'saccade peths');


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
        folderFiles = dir(fullfile(srcFull, '*.m'));
        if isempty(folderFiles)
            fprintf('  EMPTY FOLDER (no .m files): %s\n', srcRel);
            continue;
        end
        fprintf('\n--- Merging folder: %s → %s ---\n', srcRel, dstDir);
        for jj = 1:length(folderFiles)
            srcFile = fullfile(srcFull, folderFiles(jj).name);
            dstFile = fullfile(dstDir, folderFiles(jj).name);
            if exist(dstFile, 'file')
                fprintf('  CONFLICT (already exists): %s\n', folderFiles(jj).name);
                nSkipped = nSkipped + 1;
                continue;
            end
            if dryRun
                fprintf('  [DRY] %s → %s\n', folderFiles(jj).name, dstDir);
            else
                try
                    movefile(srcFile, dstFile);
                    fprintf('  MOVED: %s\n', folderFiles(jj).name);
                catch ME
                    fprintf('  FAILED: %s — %s\n', folderFiles(jj).name, ME.message);
                    nFailed = nFailed + 1;
                end
            end
            nMoved = nMoved + 1;
        end
        folderMerges{end+1} = srcFull; %#ok<SAGROW>
        continue;
    end
    
    % --- Handle deletion ---
    if strcmp(dstDir, '__DELETE__')
        if ~exist(srcFull, 'file')
            fprintf('  SKIP (not found): %s\n', srcRel);
            nSkipped = nSkipped + 1;
            continue;
        end
        if dryRun
            fprintf('[DRY] DELETE: %s\n', srcRel);
        else
            delete(srcFull);
            fprintf('DELETED: %s\n', srcRel);
        end
        nMoved = nMoved + 1;
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
    fprintf('  Would move/delete: %d files\n', nMoved);
else
    fprintf('  Moved/deleted: %d files\n', nMoved);
end
fprintf('  Skipped (not found or conflict): %d\n', nSkipped);
if ~dryRun
    fprintf('  Failed: %d\n', nFailed);
end

% --- Report empty source folders that can be deleted manually ---
if ~isempty(folderMerges)
    fprintf('\n  Folders that were merged (delete manually if empty):\n');
    for ii = 1:length(folderMerges)
        remaining = dir(fullfile(folderMerges{ii}, '*'));
        remaining = remaining(~ismember({remaining.name}, {'.', '..'}));
        if isempty(remaining)
            fprintf('    [EMPTY — safe to delete] %s\n', folderMerges{ii});
        else
            fprintf('    [%d items remain] %s\n', length(remaining), folderMerges{ii});
        end
    end
end

% --- What stays in helper functions/ after reorg ---
fprintf('\n  Files remaining in helper functions/ after reorg:\n');
hfDir = fullfile(repo, 'helper functions');
hfRemaining = dir(fullfile(hfDir, '*.m'));
% Filter out files we''re moving
movedNames = {};
for iM = 1:size(moves,1)
    [~, fn, fe] = fileparts(moves{iM,1});
    movedNames{end+1} = [fn fe]; %#ok<SAGROW>
end
for jj = 1:length(hfRemaining)
    if ~ismember(hfRemaining(jj).name, movedNames)
        fprintf('    %s\n', hfRemaining(jj).name);
    end
end

fprintf('\n');
if dryRun
    fprintf('>>> Set dryRun = false and re-run to execute. <<<\n');
else
    fprintf('>>> Done! Run CLWD_inventory_ver_1_0 to verify. <<<\n');
end
