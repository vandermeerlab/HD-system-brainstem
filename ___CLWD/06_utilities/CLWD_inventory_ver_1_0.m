function CLWD_inventory_ver_1_0()
% CLWD_inventory_ver_1_0 - List all .m files in ___CLWD and save to for_claude/
%
% Usage:
%   CLWD_inventory_ver_1_0
%
% Clears out for_claude/ folder first, then writes a text file listing
% all .m files with their subfolder locations. Output goes to:
%   ___CLWD\for_claude\CLWD_inventory.txt
%
% Run this at the start of every Claude session, then upload the .txt file.

clwd_root = 'C:\Jeff\Github_repositories\HD-system-brainstem\___CLWD';
dst = fullfile(clwd_root, 'for_claude');

% --- Create for_claude if needed, delete existing contents ---
if ~exist(dst, 'dir')
    mkdir(dst);
else
    old_files = dir(fullfile(dst, '*'));
    for iF = 1:numel(old_files)
        if ~old_files(iF).isdir
            delete(fullfile(dst, old_files(iF).name));
        end
    end
end

% --- Get all .m files recursively ---
flist = dir(fullfile(clwd_root, '**', '*.m'));

% --- Build inventory text ---
lines = cell(numel(flist) + 3, 1);
lines{1} = sprintf('=== ___CLWD inventory (%d .m files) ===', numel(flist));
lines{2} = sprintf('Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
lines{3} = '';

for ii = 1:numel(flist)
    rel_folder = strrep(flist(ii).folder, [clwd_root filesep], '');
    if strcmp(rel_folder, flist(ii).folder)
        rel_folder = '.';  % file is at ___CLWD root level
    end
    lines{ii + 3} = sprintf('  %-50s  %s', flist(ii).name, rel_folder);
end

inventory_text = strjoin(lines, newline);

% --- Print to command window ---
fprintf('\n%s\n\n', inventory_text);

% --- Save to for_claude ---
out_path = fullfile(dst, 'CLWD_inventory.txt');
fid = fopen(out_path, 'w');
fprintf(fid, '%s', inventory_text);
fclose(fid);

fprintf('Saved to: %s\n', out_path);
fprintf('Upload this file to Claude.\n');

end

% =========================================================================
% CHANGELOG (reverse chronological)
% =========================================================================
% 2026-04-04  v1.0  Initial version. Lists all .m files in ___CLWD,
%                    clears for_claude/ first, saves inventory .txt.
