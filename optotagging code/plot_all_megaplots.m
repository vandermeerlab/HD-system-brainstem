function plot_all_megaplots(tfilelist, savedir, startCell)
% JJS. 2025-09-20. Function to plot all megaplots (subplots with fingerprint for each cell) and save them in a given directory
%   Inputs
%       tfilelist: cell array where each cell is the path for a single neuron .t file
%       savedir:   path for the directory in which we want to save the files
%   Outputs
cd('C:\Jeff\U01 data\datatouse');
if isempty(tfilelist)
    tfilelist = FindFiles('*.t');
end
if isempty(savedir)
    pushdir('D:\Jeff\U01\NPH manuscript\_analysis\optotagging\megaplots');
    %     savedir = 'D:\Jeff\U01\NPH manuscript\_analysis\optotagging\megaplots';
    today = string(datetime('today'));
    mkdir(today);
    newfolder = strcat('D:\Jeff\U01\NPH manuscript\_analysis\optotagging\megaplots\', today);
    %     cd(newfolder);
    savedir = newfolder;
    popdir;
end
if isempty(startCell)
    startCell = 1;
end

for iCell = startCell:length(tfilelist)
    disp(num2str(iCell))
    [a, b, ~] = fileparts(tfilelist{iCell});  % b is a string -- name of the current neuron
    
    [~, e, ~] = fileparts(a);  % get session name for current cell
    [~, h, ~] = fileparts(pwd); % get session name for current directory
%     if strcmp(e, h) ~=1          % if sessions not the same, push to the current cell dir
%         pushdir(fileparts(tfilelist{iCell}));
%         sd = LoadSessionData([]);
%     end
    pushdir(fileparts(tfilelist{iCell}));
    sd = LoadSessionData([]);
    tlist = FindFiles('*.t');    % get list of neurons for this session
    [~, h, ~] = fileparts(tlist);
    x = strcmp(b, h);
    y = find(x);
    currentCell = y;
    
    
    cfg_in = [];
    fingerprintPlot(cfg_in, sd, currentCell)
    fig = gcf;
    fig.WindowState = 'fullscreen';
    disp('saving...')
    filename = b; % 
    fullDestination = fullfile(savedir, filename)
    saveas(gcf, fullDestination, 'png');
    clf
end

