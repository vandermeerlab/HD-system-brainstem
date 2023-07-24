function [start_time, stop_time, laser_on, laser_off, arraysize] = sessionLoopTemplate(fd, startSess, endSess)
if (isempty(fd) == 1)
    fd = FindFiles('*keys.m');
end
%     labels = [];
if isempty(startSess) ==1
    startSess = 1;
end
if isempty(endSess) ==1
    endSess = length(fd);
end
for iSess = startSess:endSess
    pushdir(fileparts(fd{iSess}));
    [start_time{iSess}, stop_time{iSess}, laser_on{iSess}, laser_off{iSess}, arraysize(iSess)] = SortBrainstemEventLabels3;
    
    popdir;
end