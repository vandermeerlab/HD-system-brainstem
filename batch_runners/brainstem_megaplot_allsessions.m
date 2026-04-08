function brainstem_megaplot_allsessions(fd, startSess, endSess)
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
%     SSN = HD_GetSSN; disp(SSN);
    sd = LoadSessionData([]);
    for iCell = 1:length(sd.S.t)
        fingerprintPlot([], sd, iCell)
        set(gcf, 'Position', get(0, 'Screensize'));
        pushdir('D:\Jeff\U01\analysis\figs\megaplots');
        WriteFig(strcat(sd.SSN, '-megaplot', num2str(iCell)))
        disp('figure saved')
        popdir;
    end
    popdir;
end
