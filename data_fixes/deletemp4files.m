% % Define the target directory
% targetDir = 'C:\Jeff\U01\datatouse';
% 
% % Enable recycling to move files to the Recycle Bin instead of permanent deletion (optional)
% % recycle('on'); 
% 
% % Get a list of all .mp4 files in the directory and all subfolders
% % The '**' pattern tells MATLAB to search subdirectories recursively
% fileList = dir(fullfile(targetDir, '**', '*VT1.Nvt'));
% 
% % Filter out any directories that might accidentally match the pattern
% fileList = fileList(~[fileList.isdir]);
% 
% % Check if any files were found
% if isempty(fileList)
%     fprintf('No .mp4 files found in %s or its subfolders.\n', targetDir);
% else
%     fprintf('Found %d files. Deleting...\n', length(fileList));
%     
%     % Loop through each file and delete it
%     for i = 1:length(fileList)
%         % Construct the full file path
%         fullFilePath = fullfile(fileList(i).folder, fileList(i).name);
%         
%         % Delete the file
%         delete(fullFilePath);
%         
%         % Optional: Print progress
%         fprintf('Deleted: %s\n', fullFilePath);
%     end
%     fprintf('Cleanup complete.\n');
% end