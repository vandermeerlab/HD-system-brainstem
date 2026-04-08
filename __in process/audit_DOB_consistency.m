function results = audit_DOB_consistency(varargin)
% 2026-03-12. JJS.
% Checks that every session for the same mouse has the same DOB (parsed
% from ExpKeys.bday in the keys file). Reports any mice with inconsistent
% or missing DOBs across sessions.
%
% Requirements: NPH_all_tfilelist.mat on path. Requires: parseSessionDate.m,
%               parseBdayFromKeysFile.m
%
% Inputs:
%           dataDir     - path to the parent data directory (default: 'C:\Jeff\U01\datatouse')
%           matfile     - filename for the tfilelist .mat file (default: 'NPH all tfilelist.mat')
%           varname     - variable name inside the .mat file (default: 'NPH_all_tfilelist')
%
% Outputs:
%           results     - struct array (one per mouse) with fields:
%                         .mouseID, .sessions, .dobs, .consistent (logical)

dataDir = 'C:\Jeff\U01\datatouse';
matfile = 'NPH all tfilelist.mat';
varname = 'NPH_all_tfilelist';
process_varargin(varargin);

%% cd to data directory
startDir = pwd;
cd(dataDir);

%% Load the neuron file list
tmp = load(matfile, varname);
tfilelist = tmp.(varname);
nNeurons = length(tfilelist);

%% Parse mouse ID and session folder from each path
mouseIDs = cell(nNeurons, 1);
sessionDirs = cell(nNeurons, 1);
sessionNames = cell(nNeurons, 1);

for iN = 1:nNeurons
    parts = strsplit(tfilelist{iN}, '\');
    sessionDirs{iN} = fullfile(parts{1:end-1});
    sessionNames{iN} = parts{end-1};
    mouseIDs{iN} = parts{end-2};
end

%% Get unique sessions
[uniqueSessions, firstIdx] = unique(sessionDirs, 'stable');
nSessions = length(uniqueSessions);

sess_mouseID = cell(nSessions, 1);
sess_sessionName = cell(nSessions, 1);
sess_dob = cell(nSessions, 1);  % DOB string from keys file

for iS = 1:nSessions
    sess_mouseID{iS} = mouseIDs{firstIdx(iS)};
    sess_sessionName{iS} = sessionNames{firstIdx(iS)};
    
    fd = uniqueSessions{iS};
    sess_dob{iS} = parseBdayFromKeysFile(fd);
    if isempty(sess_dob{iS})
        sess_dob{iS} = 'MISSING';
    end
end

%% Check consistency per mouse
uniqueMice = unique(mouseIDs, 'stable');
nMice = length(uniqueMice);

nIssues = 0;
results = struct('mouseID', {}, 'sessions', {}, 'dobs', {}, 'consistent', {});

for iM = 1:nMice
    mid = uniqueMice{iM};
    sessIdx = find(strcmp(sess_mouseID, mid));
    
    dobs = sess_dob(sessIdx);
    sessions = sess_sessionName(sessIdx);
    uniqueDobs = unique(dobs);
    
    isConsistent = (length(uniqueDobs) == 1) && ~strcmp(uniqueDobs{1}, 'MISSING');
    
    results(iM).mouseID = mid;
    results(iM).sessions = sessions;
    results(iM).dobs = dobs;
    results(iM).consistent = isConsistent;
    
    if ~isConsistent
        nIssues = nIssues + 1;
        fprintf('\n*** %s — DOB inconsistency ***\n', mid);
        for iSess = 1:length(sessIdx)
            fprintf('    %s : %s\n', sessions{iSess}, dobs{iSess});
        end
    end
end

cd(startDir);

%% Summary
fprintf('\n--- DOB Audit Summary ---\n');
fprintf('Total mice: %d\n', nMice);
fprintf('Mice with consistent DOB: %d\n', sum([results.consistent]));
fprintf('Mice with issues: %d\n', nIssues);

if nIssues == 0
    fprintf('All mice have consistent DOBs across sessions.\n');
end

end

%% Changelog
% v1.0  2026-03-12  JJS.  Initial version. Parses DOB from keys files
%                         across all sessions and checks within-mouse
%                         consistency.
