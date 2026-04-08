function bdayStr = parseBdayFromKeysFile(fd)
% 2026-03-12. JJS.
% Parses ExpKeys.bday as a date string from a keys .m file by reading the
% file as raw text. This is necessary because ExpKeys.bday = 2021-07-06
% evaluates to a number (2021 - 7 - 6 = 2008) in MATLAB, not a date.
%
% Inputs:
%           fd          - path to the session folder containing the keys file.
%                         If empty, uses the current directory.
%
% Outputs:
%           bdayStr     - date of birth as a string in 'yyyy-MM-dd' format
%                         (e.g. '2021-07-06'). Returns '' if bday cannot
%                         be found or parsed.

if nargin < 1 || isempty(fd)
    fd = pwd;
end

bdayStr = '';
origDir = pwd;
cd(fd);
try
    keysFile = FindFile('*keys.m');
    txt = fileread(keysFile);
    % Match patterns like: ExpKeys.bday = 2021-07-06  or  ExpKeys.bday = 2021-7-6
    tok = regexp(txt, 'ExpKeys\.bday\s*=\s*(\d{4})\s*-\s*(\d{1,2})\s*-\s*(\d{1,2})', 'tokens');
    if ~isempty(tok)
        yr = str2double(tok{1}{1});
        mo = str2double(tok{1}{2});
        dy = str2double(tok{1}{3});
        bdayStr = sprintf('%04d-%02d-%02d', yr, mo, dy);
    end
catch
    % could not find or parse keys file
end
cd(origDir);

end

%% Changelog
% v1.0  2026-03-12  JJS.  Initial version. Extracted from
%                         make_mouse_table_ver_1_2.
