function findSteadyPeriodsManual_AllSessions(fd, startSess, endSess)
if (isempty(fd) == 1)
    fd = FindFiles('*keys.m');
end
    if isempty(startSess) ==1
        startSess = 1;
    end
    if isempty(endSess) ==1
        endSess = length(fd);
    end
    for iSess = startSess:endSess
        pushdir(fileparts(fd{iSess}));
        [STtstart, STtend, tlistX] = findSteadyPeriodsManual2;
        popdir;
        clf
    end
end