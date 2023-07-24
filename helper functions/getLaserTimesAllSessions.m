function [START_TIME, STOP_TIME, LASER_ON, LASER_OFF, DIFFLASER, ARRAY_SIZE] = getLaserTimesAllSessions(fd, startSess, endSess)
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
    [START_TIME{iSess}, STOP_TIME{iSess}, LASER_ON{iSess}, LASER_OFF{iSess}, ARRAY_SIZE(iSess), DIFFLASER{iSess}] = SortBrainstemEventLabels3;
%     DURATION{iSESS} = LASER_OFF{iSess} - LASER_ON{iSess}; 
    popdir;
end