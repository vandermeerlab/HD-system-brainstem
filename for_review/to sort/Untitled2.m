% 1. Enable recycling so delete moves to Trash/Recycle Bin
oldState = recycle('on');

baseDir = 'C:\Jeff\U01\datatouse'; % Your directory

% 2. Search for everything inside subfolders
% Using '**/*' ensures we see every file and folder in the tree
allItems = dir(fullfile(baseDir, '**', '*'));

% 3. Filter for items whose PARENT folder is 'FD'
% This identifies everything INSIDE an 'FD' folder
isInFD = strcmp({allItems.folder}, 'FD') | ...
         cellfun(@(x) endsWith(x, [filesep 'FD']), {allItems.folder});

% 4. Collect the unique paths of the 'FD' folders themselves
% We extract the folder path for any item found inside an 'FD' directory
fdFolderPaths = unique({allItems(isInFD).folder});

% 5. Process each identified 'FD' folder
for i = 1:numel(fdFolderPaths)
    target = fdFolderPaths{i};
    
    % Verify the folder name is exactly 'FD' and it exists
    [~, name] = fileparts(target);
    if strcmp(name, 'FD') && exist(target, 'dir')
        try
            % Move contents to Recycle Bin
            % This preserves the files while 'rmdir' cleans the shell
            delete(fullfile(target, '*'));
            
            % Remove the empty 'FD' folder
            [status, msg] = rmdir(target, 's');
            
            if status
                fprintf('Recycled: %s\n', target);
            else
                fprintf('Could not remove shell %s: %s\n', target, msg);
            end
        catch ME
            fprintf('Error: %s\n', ME.message);
        end
    end
end

% Restore original recycle setting
recycle(oldState);



