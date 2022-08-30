function [nasalSaccadesToUse, temporalSaccadesToUse] = restrictSaccadesToStationaryPeriods

SSN = HD_GetSSN; disp(SSN);
load(FindFile('*AHV_StationaryTimes.mat'), 'STtstart', 'STtend'); 
load(FindFile(strcat(SSN, '-saccades-edited', '.mat')));


N.t = {nasalSaccades};
T.t = {temporalSaccades};

nasalSaccadesToUse = restrict(N, STtstart, STtend);
temporalSaccadesToUse = restrict(T, STtstart, STtend);

