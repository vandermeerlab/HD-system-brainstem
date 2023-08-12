function [] = save_a_copy(fd)

if isempty(fd)
    fd = FindFiles('*keys.m');
end
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
%     SSN = HD_GetSSN; disp(SSN);
    [START_TIME{iSess}, STOP_TIME{iSess}, LASER_ON{iSess}, LASER_OFF{iSess}, BIT0{iSess}, BIT4{iSess}, ARRAYSIZE(iSess)] = SortBrainstemEventLabels2;
end
