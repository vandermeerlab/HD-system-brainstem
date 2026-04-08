function sessDate = parseSessionDate(sessName)
% 2026-03-12. JJS.
% Parses the recording date from a session folder name.
%
% Inputs:
%           sessName    - session folder name string 
%                         e.g. 'M039-2020-08-21-1' or 'M444-2023-08-18'
%
% Outputs:
%           sessDate    - datetime object for the session date.
%                         Returns NaT if the date cannot be parsed.

tok = regexp(sessName, 'M\d+-(\d{4}-\d{2}-\d{2})', 'tokens');
if ~isempty(tok)
    sessDate = datetime(tok{1}{1}, 'InputFormat', 'yyyy-MM-dd');
else
    sessDate = NaT;
end

end

%% Changelog
% v1.0  2026-03-12  JJS.  Initial version. Extracted from
%                         make_mouse_table_ver_1_2.
