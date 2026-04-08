function NPH_organize_claude_dir_ver_1_0()
% NPH_organize_claude_dir_ver_1_0
%
% Organizes the ___CLWD folder into categorical subfolders and moves older
% versioned files to the 'old_versions' subfolder.
%
% Run this at the start of each session to keep the folder tidy.
% Files are MOVED (not copied) so there are no duplicates.
%
% CATEGORIES (subfolders):
%   01_data_loading/     - LoadSpikesQuick
%   02_EV_pipeline/      - compute_slow_phase_EV, diagnose_desaccade_order
%   03_screening/        - screen_AHV_vs_EV_opposition
%   04_optotag/          - compute_optotag_shuffle, run_optotag_shuffle,
%                          plot_optotag_laser_PETHs
%   05_figure_code/      - plot_fig4B_*, plot_AHV_vs_EV_figure,
%                          plot_EV_TC_from_workspace, plot_EV_gappy_vs_interp
%   06_utilities/        - generate_all_function_paths, start_of_session_copy,
%                          NPH_organize_claude_dir
%   old_versions/        - Superseded versions
%   for_claude/          - Session staging area (untouched)
%
% USAGE:
%   NPH_organize_claude_dir_ver_1_0     % runs with confirmation prompt
%
% AFTER RUNNING:
%   restoredefaultpath; clear classes;
%   addpath(genpath('C:\Jeff\Github_repositories\vandermeerlab\code-matlab\shared'))
%   addpath(genpath('C:\Jeff\Github_repositories\HD-system-brainstem'))
%
% Version: 1.0
% 2026-04-03

%% ========== CONFIGURATION ==========
CLAUDE_DIR = 'C:\Jeff\Github_repositories\HD-system-brainstem\___CLWD';
OLD_DIR_NAME = 'old_versions';

%% ========== DEFINE CATEGORY ASSIGNMENTS ==========
% Each row: { 'filename.m', 'subfolder_name' }
% Only the LATEST version of each function goes here.
% Files not listed will be flagged as unassigned and LEFT IN PLACE.

assignments = {
    % --- 01_data_loading ---
    'LoadSpikesQuick.m',                        '01_data_loading'

    % --- 02_EV_pipeline ---
    'compute_slow_phase_EV.m',                  '02_EV_pipeline'
    'diagnose_desaccade_order.m',               '02_EV_pipeline'

    % --- 03_screening ---
    'screen_AHV_vs_EV_opposition.m',            '03_screening'

    % --- 04_optotag ---
    'compute_optotag_shuffle.m',                '04_optotag'
    'run_optotag_shuffle.m',                    '04_optotag'
    'plot_optotag_laser_PETHs.m',               '04_optotag'

    % --- 05_figure_code ---
    'plot_fig4B_ver_2_1.m',                     '05_figure_code'
    'plot_AHV_vs_EV_figure.m',                  '05_figure_code'
    'plot_EV_TC_from_workspace.m',              '05_figure_code'
    'plot_EV_gappy_vs_interp.m',                '05_figure_code'

    % --- 06_utilities ---
    'generate_all_function_paths.m',            '06_utilities'
    'start_of_session_copy.m',                  '06_utilities'
    'NPH_organize_claude_dir_ver_1_0.m',        '06_utilities'
};

%% ========== AUTO-DETECT OLD VERSIONS ==========
% Scan top-level for versioned files and identify superseded ones.
% A file is "old" if a newer version of the same base name exists.

topLevelListing = dir(fullfile(CLAUDE_DIR, '*.m'));
topLevelNames = {topLevelListing.name}';
assignedFiles = assignments(:,1);

% Parse versioned filenames: baseName_ver_X_Y.m
oldVersionFiles = {};
for ii = 1:length(topLevelNames)
    fname = topLevelNames{ii};
    % Check if this file is already assigned to a category
    if any(strcmp(fname, assignedFiles))
        continue  % it's a current file, skip
    end
    % Check if it matches a versioned pattern of an assigned file
    % e.g., plot_fig4B_ver_1_5.m is an old version of plot_fig4B_ver_2_1.m
    tok = regexp(fname, '^(.+)_ver_(\d+)_(\d+)\.m$', 'tokens');
    if ~isempty(tok)
        baseName = tok{1}{1};
        % See if any assigned file shares this base name
        for jj = 1:length(assignedFiles)
            tok2 = regexp(assignedFiles{jj}, '^(.+)_ver_(\d+)_(\d+)\.m$', 'tokens');
            if ~isempty(tok2) && strcmp(tok2{1}{1}, baseName)
                oldVersionFiles{end+1} = fname; %#ok<AGROW>
                break
            end
        end
    end
end
oldVersionFiles = oldVersionFiles(:);

%% ========== VALIDATE CLAUDE_DIR ==========
if ~isfolder(CLAUDE_DIR)
    error('CLAUDE_DIR not found: %s', CLAUDE_DIR);
end

oldDir = fullfile(CLAUDE_DIR, OLD_DIR_NAME);
if ~isfolder(oldDir)
    mkdir(oldDir);
    fprintf('Created: %s\n', OLD_DIR_NAME);
end

assignedFolders = assignments(:,2);
subfolderNames = unique(assignedFolders);

%% ========== DETECT UNASSIGNED FILES ==========
allKnownFiles = [assignedFiles; oldVersionFiles];
unassignedFiles = setdiff(topLevelNames, allKnownFiles);

% Filter out process_varargin (shared dependency, intentionally at top level)
keepAtTopLevel = {'process_varargin.m'};
unassignedFiles = setdiff(unassignedFiles, keepAtTopLevel);

if ~isempty(unassignedFiles)
    fprintf('NOTE: %d file(s) not assigned to any category:\n', length(unassignedFiles));
    for ii = 1:length(unassignedFiles)
        fprintf('  ?? %s\n', unassignedFiles{ii});
    end
    fprintf('  These files will be LEFT IN PLACE at the top level.\n\n');
end

%% ========== COUNT WHAT NEEDS MOVING ==========
% Category moves: files at top level that should be in subfolders
nMoveCategory = 0;
for ii = 1:size(assignments, 1)
    srcPath = fullfile(CLAUDE_DIR, assignments{ii, 1});
    dstPath = fullfile(CLAUDE_DIR, assignments{ii, 2}, assignments{ii, 1});
    if isfile(srcPath) && ~isfile(dstPath)
        nMoveCategory = nMoveCategory + 1;
    end
end

% Old version moves
oldExist = {};
for ii = 1:length(oldVersionFiles)
    srcPath = fullfile(CLAUDE_DIR, oldVersionFiles{ii});
    if isfile(srcPath)
        oldExist{end+1} = oldVersionFiles{ii}; %#ok<AGROW>
    end
end
nMoveOld = length(oldExist);

totalMoves = nMoveCategory + nMoveOld;

if totalMoves == 0
    fprintf('___CLWD is already organized. Nothing to move.\n');
    return
end

%% ========== PRINT PLAN ==========
fprintf('============================================================\n');
fprintf('  ___CLWD ORGANIZATION PLAN\n');
fprintf('============================================================\n\n');

for iFolder = 1:length(subfolderNames)
    folderName = subfolderNames{iFolder};
    folderIdx = strcmp(assignedFolders, folderName);
    filesInFolder = assignedFiles(folderIdx);
    filesExist = {};
    for jj = 1:length(filesInFolder)
        srcPath = fullfile(CLAUDE_DIR, filesInFolder{jj});
        dstPath = fullfile(CLAUDE_DIR, folderName, filesInFolder{jj});
        if isfile(srcPath) && ~isfile(dstPath)
            filesExist{end+1} = filesInFolder{jj}; %#ok<AGROW>
        end
    end
    if isempty(filesExist); continue; end
    fprintf('  %s/  (%d files)\n', folderName, length(filesExist));
    for jj = 1:length(filesExist)
        fprintf('    -> %s\n', filesExist{jj});
    end
    fprintf('\n');
end

if nMoveOld > 0
    fprintf('  %s/  (%d old versions)\n', OLD_DIR_NAME, nMoveOld);
    for jj = 1:nMoveOld
        fprintf('    -> %s\n', oldExist{jj});
    end
    fprintf('\n');
end

fprintf('------------------------------------------------------------\n');
fprintf('  TOTAL: %d files to move\n', totalMoves);
fprintf('------------------------------------------------------------\n\n');

%% ========== CONFIRM ==========
confirmStr = input('Type YES to proceed: ', 's');
if ~strcmpi(strtrim(confirmStr), 'YES')
    fprintf('Aborted.\n');
    return
end

%% ========== CREATE SUBFOLDERS ==========
for iFolder = 1:length(subfolderNames)
    folderPath = fullfile(CLAUDE_DIR, subfolderNames{iFolder});
    if ~isfolder(folderPath)
        mkdir(folderPath);
        fprintf('Created: %s\n', subfolderNames{iFolder});
    end
end

%% ========== MOVE OLD VERSIONS ==========
nMovedOld = 0;
for ii = 1:nMoveOld
    srcPath = fullfile(CLAUDE_DIR, oldExist{ii});
    dstPath = fullfile(oldDir, oldExist{ii});
    if isfile(dstPath)
        fprintf('  SKIP (already in %s): %s\n', OLD_DIR_NAME, oldExist{ii});
        continue
    end
    [moveOK, moveMsg] = movefile(srcPath, dstPath);
    if moveOK
        fprintf('  -> %s: %s\n', OLD_DIR_NAME, oldExist{ii});
        nMovedOld = nMovedOld + 1;
    else
        fprintf('  FAILED: %s -- %s\n', oldExist{ii}, moveMsg);
    end
end

%% ========== MOVE CURRENT FILES TO SUBFOLDERS ==========
nMovedCategory = 0;
for ii = 1:size(assignments, 1)
    fileName = assignments{ii, 1};
    folderName = assignments{ii, 2};
    srcPath = fullfile(CLAUDE_DIR, fileName);
    dstPath = fullfile(CLAUDE_DIR, folderName, fileName);

    if ~isfile(srcPath)
        continue
    end
    if isfile(dstPath)
        % File exists in both places — remove top-level duplicate
        delete(srcPath);
        fprintf('  Removed duplicate: %s (already in %s)\n', fileName, folderName);
        nMovedCategory = nMovedCategory + 1;
        continue
    end
    [moveOK, moveMsg] = movefile(srcPath, dstPath);
    if moveOK
        fprintf('  -> %s: %s\n', folderName, fileName);
        nMovedCategory = nMovedCategory + 1;
    else
        fprintf('  FAILED: %s -- %s\n', fileName, moveMsg);
    end
end

%% ========== FINAL REPORT ==========
fprintf('\n============================================================\n');
fprintf('  DONE\n');
fprintf('  %d files -> category subfolders\n', nMovedCategory);
fprintf('  %d files -> %s\n', nMovedOld, OLD_DIR_NAME);
fprintf('============================================================\n\n');
fprintf('  Refresh your MATLAB path:\n\n');
fprintf('    restoredefaultpath; clear classes;\n');
fprintf('    addpath(genpath(''C:\\Jeff\\Github_repositories\\vandermeerlab\\code-matlab\\shared''))\n');
fprintf('    addpath(genpath(''C:\\Jeff\\Github_repositories\\HD-system-brainstem''))\n\n');

end

%% Version History
% 2026-04-03  v1.0  Initial version. Adapted from PM_organize_claude_dir_ver_1_0
%                    for the NPH/brainstem project ___CLWD folder.
