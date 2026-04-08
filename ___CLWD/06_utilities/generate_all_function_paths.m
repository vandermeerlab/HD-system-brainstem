% generate_all_function_paths.m
% Save this as a .m file and run it. Do not paste into command window.
%
% Scans all 4 repos on Jeff's MATLAB path and writes every .m file
% with its full path to a text file for Claude's project KB.

repos = {
    'C:\Jeff\Github_repositories\vandermeerlab\code-matlab\shared', 'vandermeerlab-shared';
    'C:\Jeff\Github_repositories\HD-system-brainstem', 'HD-system-brainstem';
    'C:\Jeff\Github_repositories\matlas\Jefflas\brainstem_registered_coordinates', 'matlas-brainstem';
    'C:\Jeff\Github_repositories\vandermeerlab\code-matlab\tasks\Alyssa_Tmaze', 'Alyssa-Tmaze';
};

out_file = 'C:\Jeff\Github_repositories\HD-system-brainstem\___CLWD\for_claude\all_function_paths.txt';

all_lines = {};
all_lines{end+1} = 'Master Function Path List';
all_lines{end+1} = ['Generated: ' datestr(now, 'yyyy-mm-dd')];
all_lines{end+1} = '';

total_count = 0;
for iR = 1:size(repos, 1)
    repo_path = repos{iR, 1};
    repo_label = repos{iR, 2};
    flist = dir(fullfile(repo_path, '**', '*.m'));
    all_lines{end+1} = ['=== ' repo_label ' (' num2str(numel(flist)) ' files) ===']; %#ok<*SAGROW>
    all_lines{end+1} = ['Root: ' repo_path];
    all_lines{end+1} = '';
    for iF = 1:numel(flist)
        all_lines{end+1} = fullfile(flist(iF).folder, flist(iF).name);
    end
    all_lines{end+1} = '';
    total_count = total_count + numel(flist);
    disp([repo_label ': ' num2str(numel(flist)) ' files']);
end

% Write all at once
fid = fopen(out_file, 'w');
for iL = 1:numel(all_lines)
    fwrite(fid, [all_lines{iL} newline]);
end
fclose(fid);

disp(' ');
disp(['Total: ' num2str(total_count) ' .m files']);
disp(['Written to: ' out_file]);
