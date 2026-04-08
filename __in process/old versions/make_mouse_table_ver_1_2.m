function mouseTable = make_mouse_table_ver_1_2(varargin)
% 2026-03-12. JJS.
% Builds a summary table of mouse metadata (sex, genotype/strain, age at
% recording, number of sessions, number of neurons) for the primary NPH
% dataset (n = 104 neurons).
%
% Requirements: NPH_all_tfilelist.mat must be on the path. Must be able to
%               cd into each session folder and run EvalKeys.
%               Requires: loadKeysFromDir.m, parseSessionDate.m,
%               parseBdayFromKeysFile.m
%
% Inputs:
%           dataDir     - path to the parent data directory (default: 'C:\Jeff\U01\datatouse')
%           matfile     - filename for the tfilelist .mat file (default: 'NPH all tfilelist.mat')
%           varname     - variable name inside the .mat file (default: 'NPH_all_tfilelist')
%           ageTolerance - max allowed discrepancy (days) between ExpKeys.age and
%                          calculated age before throwing a warning (default: 2)
%
% Outputs:
%           mouseTable  - MATLAB table with columns: MouseID, Sex, Strain,
%                         AgeMin_days, AgeMax_days, NumSessions, NumNeurons

dataDir = 'C:\Jeff\U01\datatouse';
matfile = 'NPH all tfilelist.mat';
varname = 'NPH_all_tfilelist';
ageTolerance = 2;
process_varargin(varargin);

%% cd to data directory
startDir = pwd;
cd(dataDir);

%% Load the neuron file list
tmp = load(matfile, varname);
tfilelist = tmp.(varname);
nNeurons = length(tfilelist);
fprintf('Total neurons in list: %d\n', nNeurons);

%% Parse mouse ID and session folder from each path
mouseIDs = cell(nNeurons, 1);
sessionDirs = cell(nNeurons, 1);
sessionNames = cell(nNeurons, 1);

for iN = 1:nNeurons
    parts = strsplit(tfilelist{iN}, '\');
    sessionDirs{iN} = fullfile(parts{1:end-1});   % full path to session folder
    sessionNames{iN} = parts{end-1};               % e.g. 'M039-2020-08-21-1'
    mouseIDs{iN} = parts{end-2};                    % e.g. 'M039'
end

%% Get unique sessions and loop through to collect metadata
[uniqueSessions, firstIdx] = unique(sessionDirs, 'stable');
nSessions = length(uniqueSessions);
fprintf('Total unique sessions: %d\n', nSessions);

% Preallocate per-session metadata
sess_mouseID = cell(nSessions, 1);
sess_sessionName = cell(nSessions, 1);
sess_sex = cell(nSessions, 1);
sess_strain = cell(nSessions, 1);
sess_age = nan(nSessions, 1);
sess_nNeurons = zeros(nSessions, 1);

for iS = 1:nSessions
    fd = uniqueSessions{iS};
    sess_sessionName{iS} = sessionNames{firstIdx(iS)};
    
    % Count neurons in this session
    sess_nNeurons(iS) = sum(strcmp(sessionDirs, fd));
    sess_mouseID{iS} = mouseIDs{firstIdx(iS)};
    
    % Parse session date from folder name
    sessDate = parseSessionDate(sess_sessionName{iS});
    
    try
        keys = loadKeysFromDir(fd);
        
        if isempty(keys.sex); sess_sex{iS} = 'unknown'; else; sess_sex{iS} = keys.sex; end
        if isempty(keys.strain); sess_strain{iS} = 'unknown'; else; sess_strain{iS} = keys.strain; end
        
        % --- Age calculation and cross-check ---
        calcAge = NaN;
        bdayStr = parseBdayFromKeysFile(fd);
        if ~isempty(bdayStr) && ~isnat(sessDate)
            try
                dob = datetime(bdayStr, 'InputFormat', 'yyyy-MM-dd');
                calcAge = days(sessDate - dob);
            catch
                % could not parse DOB string
            end
        end
        
        keysAge = [];
        if isfield(keys, 'age'); keysAge = keys.age; end
        
        if ~isempty(keysAge) && ~isnan(calcAge)
            % Both exist: cross-check
            if abs(keysAge - calcAge) > ageTolerance
                warning('AGE MISMATCH: %s — ExpKeys.age = %d, calculated from DOB = %d (diff = %d days)', ...
                    sess_sessionName{iS}, keysAge, calcAge, keysAge - calcAge);
            end
            sess_age(iS) = keysAge;  % trust ExpKeys.age when both are present
        elseif ~isempty(keysAge)
            % Only ExpKeys.age exists
            sess_age(iS) = keysAge;
        elseif ~isnan(calcAge)
            % Only calculated age exists — fill in the gap
            fprintf('  Age calculated from DOB for %s: %d days\n', sess_sessionName{iS}, calcAge);
            sess_age(iS) = calcAge;
        else
            sess_age(iS) = NaN;
        end
        
    catch ME
        warning('Could not load keys for session %s: %s', fd, ME.message);
        sess_sex{iS} = 'unknown';
        sess_strain{iS} = 'unknown';
    end
end
cd(startDir);

%% Normalize strain names
sess_strain = strrep(sess_strain, 'Ai-32', 'Ai32');

%% Aggregate per mouse
uniqueMice = unique(mouseIDs, 'stable');
nMice = length(uniqueMice);

mouse_sex = cell(nMice, 1);
mouse_strain = cell(nMice, 1);
mouse_minAge = nan(nMice, 1);
mouse_maxAge = nan(nMice, 1);
mouse_nSessions = zeros(nMice, 1);
mouse_nNeurons = zeros(nMice, 1);

for iM = 1:nMice
    mid = uniqueMice{iM};
    
    % Find sessions belonging to this mouse
    sessIdx = find(strcmp(sess_mouseID, mid));
    mouse_nSessions(iM) = length(sessIdx);
    mouse_nNeurons(iM) = sum(sess_nNeurons(sessIdx));
    
    % Sex and strain should be consistent within a mouse — take the first
    mouse_sex{iM} = sess_sex{sessIdx(1)};
    mouse_strain{iM} = sess_strain{sessIdx(1)};
    
    % Age range across recording sessions
    ages = sess_age(sessIdx);
    mouse_minAge(iM) = min(ages);
    mouse_maxAge(iM) = max(ages);
end

%% Build the table
mouseTable = table(uniqueMice, mouse_sex, mouse_strain, mouse_minAge, mouse_maxAge, mouse_nSessions, mouse_nNeurons, ...
    'VariableNames', {'MouseID', 'Sex', 'Strain', 'AgeMin_days', 'AgeMax_days', 'NumSessions', 'NumNeurons'});

%% Display
disp(mouseTable);
fprintf('\n--- Summary ---\n');
fprintf('Total mice: %d\n', nMice);
fprintf('Total sessions: %d\n', sum(mouse_nSessions));
fprintf('Total neurons: %d\n', sum(mouse_nNeurons));

% Sex breakdown
nMale = sum(strcmpi(mouse_sex, 'Male'));
nFemale = sum(strcmpi(mouse_sex, 'Female'));
fprintf('Males: %d, Females: %d\n', nMale, nFemale);

% Strain breakdown
strains = unique(mouse_strain);
for iStr = 1:length(strains)
    n = sum(strcmp(mouse_strain, strains{iStr}));
    fprintf('Strain %s: %d mice\n', strains{iStr}, n);
end

% Age range across all mice
fprintf('Age range: %d to %d days (%.1f to %.1f months)\n', ...
    min(mouse_minAge), max(mouse_maxAge), min(mouse_minAge)/30.4, max(mouse_maxAge)/30.4);

end

%% Changelog
% v1.2  2026-03-12  JJS.  Added automatic age calculation from DOB
%                         (parsed as text from keys file via
%                         parseBdayFromKeysFile.m) and session date
%                         (parsed from folder name via parseSessionDate.m).
%                         If ExpKeys.age is empty but DOB is present, age
%                         is filled in automatically. If both exist,
%                         cross-checks and warns if discrepancy exceeds
%                         ageTolerance (default 2 days). Extracted all
%                         helper subfunctions to standalone files:
%                         loadKeysFromDir.m, parseSessionDate.m,
%                         parseBdayFromKeysFile.m.
% v1.1  2026-03-12  JJS.  Moved clear ExpKeys to before EvalKeys call to
%                         prevent stale struct fields from causing size
%                         mismatch errors (fixed M284, M287, M293).
%                         Added strain normalization ('Ai-32' -> 'Ai32').
%                         Added dataDir parameter (default: C:\Jeff\U01\datatouse).
%                         Converted to function with process_varargin.
%                         Added changelog.
% v1.0  2026-03-12  JJS.  Initial version (script). Loops through
%                         NPH_all_tfilelist, loads ExpKeys via EvalKeys,
%                         extracts sex/strain/age, aggregates per mouse.
