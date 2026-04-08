function keys = loadKeysFromDir(fd)
% 2026-03-12. JJS.
% Loads ExpKeys from a session directory in an isolated function workspace.
% This prevents stale ExpKeys fields from a previous session from causing
% size mismatch errors when EvalKeys overwrites struct fields.
%
% Inputs:
%           fd          - path to the session folder containing the keys file.
%                         If empty, uses the current directory.
%
% Outputs:
%           keys        - the ExpKeys struct loaded from the keys file.

if nargin < 1 || isempty(fd)
    fd = pwd;
end

origDir = pwd;
cd(fd);
EvalKeys;  % creates ExpKeys in this function's workspace — no conflict with caller
keys = ExpKeys;
cd(origDir);

end

%% Changelog
% v1.0  2026-03-12  JJS.  Initial version. Extracted from
%                         make_mouse_table_ver_1_2.
