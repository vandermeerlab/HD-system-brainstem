function [C, index, sessList] = CountNumNeuronsPerSession(fd, minCells)
% fd = the directory of sessions 
% minCells = min number of neurons for a session that you'd like to select for
%       if the # of neurons is unimportant, set threshold = 0;

% fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    c = FindFiles('*.t');
    num = length(c);
    C(iSess) = num;
    popdir;
end

[index] = find(C >= minCells);
sessList = fd{index};

% SOME weird bug here. fd{index} is not working, and it should 

% if isempty(fd)
%     fd = FindFiles('*keys.m');
%     for iSess = 1:length(fd)
%         fd{iSess} = fileparts(fd{iSess});
%     end
% end





