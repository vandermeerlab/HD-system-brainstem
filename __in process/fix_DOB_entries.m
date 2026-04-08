%% fix_DOB_entries.m
% 2026-03-12. JJS.
% One-off script to correct ExpKeys.bday entries in keys files.
% Reads each keys file as text, finds/replaces or inserts the bday line,
% and writes the corrected file back to disk.
%
% Run once from any directory. Backs up original keys file as *_keys_backup.m
% before overwriting.

dataDir = 'C:\Jeff\U01\datatouse';

% Define fixes: {session folder, correct DOB string}
% For sessions with an existing bday line, it will be replaced.
% For sessions missing bday, the line will be inserted after ExpKeys.age.
fixes = {
    'M039\M039-2020-08-21-2',   '2019-11-13'   % was 2019-11-03
    'M039\M039-2020-08-22',     '2019-11-13'   % was 2019-11-03
    'M281\M281-2021-12-23',     '2021-07-07'   % was 2021-07-21
    'M281\M281-2021-12-27',     '2021-07-07'   % was 2021-07-03
    'M281\M281-2021-12-28',     '2021-07-07'   % was 2021-07-03
    'M284\M284-2022-05-11-2',   '2021-04-27'   % was MISSING
    'M287\M287-2022-08-05',     '2021-07-06'   % was MISSING
    'M293\M293-2021-08-13',     '2021-05-05'   % was MISSING
    'M293\M293-2021-08-19',     '2021-05-05'   % was MISSING
    };

nFixes = size(fixes, 1);
fprintf('Applying %d DOB fixes...\n', nFixes);

for iF = 1:nFixes
    sessPath = fullfile(dataDir, fixes{iF, 1});
    correctDOB = fixes{iF, 2};
    
    % Find the keys file
    origDir = pwd;
    cd(sessPath);
    keysFile = FindFile('*keys.m');
    cd(origDir);
    
    % Read the file
    txt = fileread(keysFile);
    
    % Backup original
    [fdir, fname, fext] = fileparts(keysFile);
    backupFile = fullfile(fdir, [fname '_backup' fext]);
    if ~exist(backupFile, 'file')
        copyfile(keysFile, backupFile);
        fprintf('  Backed up: %s\n', backupFile);
    end
    
    % Parse the correct DOB into the YYYY-MM-DD format used in keys files
    newBdayLine = sprintf('ExpKeys.bday = %s;', correctDOB);
    
    % Check if a bday line already exists
    bdayPattern = 'ExpKeys\.bday\s*=\s*[^;]+;[^\n]*';
    if ~isempty(regexp(txt, bdayPattern, 'once'))
        % Replace existing bday line
        oldMatch = regexp(txt, bdayPattern, 'match', 'once');
        txt = regexprep(txt, bdayPattern, newBdayLine);
        fprintf('  %s: replaced "%s" with "%s"\n', fixes{iF, 1}, strtrim(oldMatch), newBdayLine);
    else
        % No bday line — insert after ExpKeys.age line, or at the end if age is also missing
        agePattern = '(ExpKeys\.age\s*=\s*[^;]+;[^\n]*)';
        if ~isempty(regexp(txt, agePattern, 'once'))
            txt = regexprep(txt, agePattern, ['$1\n' newBdayLine ' % mouse date of birth  % YYYY-MM-DD']);
            fprintf('  %s: inserted "%s" after ExpKeys.age\n', fixes{iF, 1}, newBdayLine);
        else
            % Append to end of file
            txt = [txt, sprintf('\n%s %% mouse date of birth  %% YYYY-MM-DD\n', newBdayLine)];
            fprintf('  %s: appended "%s" to end of file\n', fixes{iF, 1}, newBdayLine);
        end
    end
    
    % Write corrected file
    fid = fopen(keysFile, 'w');
    fprintf(fid, '%s', txt);
    fclose(fid);
end

fprintf('\nDone. Re-run audit_DOB_consistency to verify.\n');
