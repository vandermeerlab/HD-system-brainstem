function results = fix_ages_from_DOB(varargin)
% 2026-03-12. JJS.
% Recalculates ExpKeys.age from ExpKeys.bday (DOB) and the session date
% (parsed from the folder name). If the existing ExpKeys.age disagrees with
% the calculated age, or if ExpKeys.age is empty, the keys file is corrected
% on disk and a note is added to ExpKeys.notes.
%
% DOBs are treated as ground truth. Run audit_DOB_consistency first to
% ensure all DOBs are correct and consistent across sessions.
%
% Requirements: NPH_all_tfilelist.mat on path. Requires: parseSessionDate.m,
%               parseBdayFromKeysFile.m
%
% Inputs:
%           dataDir     - path to the parent data directory (default: 'C:\Jeff\U01\datatouse')
%           matfile     - filename for the tfilelist .mat file (default: 'NPH all tfilelist.mat')
%           varname     - variable name inside the .mat file (default: 'NPH_all_tfilelist')
%           dryRun      - if 1, report what would change but don't write files (default: 0)
%
% Outputs:
%           results     - struct array with fields: .session, .oldAge,
%                         .newAge, .action ('corrected', 'inserted', 'ok', 'no_dob')

dataDir = 'C:\Jeff\U01\datatouse';
matfile = 'NPH all tfilelist.mat';
varname = 'NPH_all_tfilelist';
dryRun = 0;
process_varargin(varargin);

%% cd to data directory
startDir = pwd;
cd(dataDir);

%% Load the neuron file list
tmp = load(matfile, varname);
tfilelist = tmp.(varname);
nNeurons = length(tfilelist);

%% Parse session folders
sessionDirs = cell(nNeurons, 1);
sessionNames = cell(nNeurons, 1);
for iN = 1:nNeurons
    parts = strsplit(tfilelist{iN}, '\');
    sessionDirs{iN} = fullfile(parts{1:end-1});
    sessionNames{iN} = parts{end-1};
end

[uniqueSessions, firstIdx] = unique(sessionDirs, 'stable');
nSessions = length(uniqueSessions);

fprintf('Checking ages for %d sessions...\n', nSessions);
if dryRun; fprintf('*** DRY RUN — no files will be modified ***\n'); end

results = struct('session', {}, 'oldAge', {}, 'newAge', {}, 'action', {});
nCorrected = 0;
nInserted = 0;
nOK = 0;
nNoDOB = 0;

for iS = 1:nSessions
    fd = uniqueSessions{iS};
    sessName = sessionNames{firstIdx(iS)};
    
    % Parse session date from folder name
    sessDate = parseSessionDate(sessName);
    if isnat(sessDate)
        warning('Could not parse date from session name: %s', sessName);
        results(iS).session = sessName;
        results(iS).oldAge = NaN;
        results(iS).newAge = NaN;
        results(iS).action = 'parse_error';
        continue
    end
    
    % Parse DOB from keys file
    bdayStr = parseBdayFromKeysFile(fd);
    if isempty(bdayStr)
        results(iS).session = sessName;
        results(iS).oldAge = NaN;
        results(iS).newAge = NaN;
        results(iS).action = 'no_dob';
        nNoDOB = nNoDOB + 1;
        continue
    end
    
    dob = datetime(bdayStr, 'InputFormat', 'yyyy-MM-dd');
    calcAge = days(sessDate - dob);
    
    % Read the keys file
    origDir = pwd;
    cd(fd);
    keysFile = FindFile('*keys.m');
    cd(origDir);
    txt = fileread(keysFile);
    
    % Find the existing ExpKeys.age value
    ageTok = regexp(txt, 'ExpKeys\.age\s*=\s*([^;]+);', 'tokens');
    if ~isempty(ageTok)
        oldAgeStr = strtrim(ageTok{1}{1});
        if strcmp(oldAgeStr, '[]')
            oldAge = [];
        else
            oldAge = str2double(oldAgeStr);
        end
    else
        oldAge = [];
    end
    
    results(iS).session = sessName;
    results(iS).oldAge = oldAge;
    results(iS).newAge = calcAge;
    
    % Decide what to do
    if ~isempty(oldAge) && oldAge == calcAge
        results(iS).action = 'ok';
        nOK = nOK + 1;
        continue
    end
    
    % Age needs fixing
    if isempty(oldAge)
        actionStr = 'inserted';
        nInserted = nInserted + 1;
        fprintf('  %s: ExpKeys.age was empty, inserting %d\n', sessName, calcAge);
    else
        actionStr = 'corrected';
        nCorrected = nCorrected + 1;
        fprintf('  %s: ExpKeys.age was %d, correcting to %d\n', sessName, oldAge, calcAge);
    end
    results(iS).action = actionStr;
    
    if dryRun; continue; end
    
    % --- Modify the keys file ---
    
    % Backup if not already done
    [fdir, fname, fext] = fileparts(keysFile);
    backupFile = fullfile(fdir, [fname '_backup' fext]);
    if ~exist(backupFile, 'file')
        copyfile(keysFile, backupFile);
    end
    
    % Replace or insert the age line
    newAgeLine = sprintf('ExpKeys.age = %d; %% in days', calcAge);
    agePattern = 'ExpKeys\.age\s*=\s*[^;]+;[^\n]*';
    if ~isempty(regexp(txt, agePattern, 'once'))
        txt = regexprep(txt, agePattern, newAgeLine);
    else
        % No age line exists — insert after bday line
        bdayPattern = '(ExpKeys\.bday\s*=\s*[^;]+;[^\n]*)';
        if ~isempty(regexp(txt, bdayPattern, 'once'))
            txt = regexprep(txt, bdayPattern, ['$1\n' newAgeLine]);
        else
            txt = [txt, sprintf('\n%s\n', newAgeLine)];
        end
    end
    
    % Add a note about the correction
    dateStamp = datestr(now, 'yyyy-mm-dd');
    noteText = sprintf('Age corrected from %s to %d (calculated from DOB %s and session date %s) on %s.', ...
        mat2str(oldAge), calcAge, bdayStr, datestr(sessDate, 'yyyy-mm-dd'), dateStamp);
    
    notesPattern = 'ExpKeys\.notes\s*=\s*\{([^}]*)\}';
    notesTok = regexp(txt, notesPattern, 'tokens');
    if ~isempty(notesTok)
        existingNotes = strtrim(notesTok{1}{1});
        if isempty(existingNotes)
            newNotesContent = sprintf('''%s''', noteText);
        else
            newNotesContent = sprintf('%s; ''%s''', existingNotes, noteText);
        end
        txt = regexprep(txt, notesPattern, sprintf('ExpKeys.notes = {%s}', newNotesContent));
    else
        % No notes field — append one
        txt = [txt, sprintf('\nExpKeys.notes = {''%s''};\n', noteText)];
    end
    
    % Write corrected file
    fid = fopen(keysFile, 'w');
    fprintf(fid, '%s', txt);
    fclose(fid);
end

cd(startDir);

%% Summary
fprintf('\n--- Age Correction Summary ---\n');
fprintf('Total sessions: %d\n', nSessions);
fprintf('Already correct: %d\n', nOK);
fprintf('Corrected: %d\n', nCorrected);
fprintf('Inserted (was empty): %d\n', nInserted);
fprintf('No DOB available: %d\n', nNoDOB);
if dryRun; fprintf('*** This was a dry run — no files were modified ***\n'); end

end

%% Changelog
% v1.0  2026-03-12  JJS.  Initial version. Recalculates ExpKeys.age from
%                         DOB and session date. Corrects keys files on disk
%                         if discrepant. Adds correction note to
%                         ExpKeys.notes. Supports dryRun mode.
