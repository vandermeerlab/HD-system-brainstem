function [START_TIME, STOP_TIME, LASER_ON, LASER_OFF, BIT0, BIT4, ARRAYSIZE] = run_laserEvents_all_sessions(fd)

if isempty(fd)
    fd = FindFiles('*keys.m');
end
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
%     SSN = HD_GetSSN; disp(SSN);
    [START_TIME{iSess}, STOP_TIME{iSess}, LASER_ON{iSess}, LASER_OFF{iSess}, BIT0{iSess}, BIT4{iSess}, ARRAYSIZE(iSess)] = SortBrainstemEventLabels2;
end
