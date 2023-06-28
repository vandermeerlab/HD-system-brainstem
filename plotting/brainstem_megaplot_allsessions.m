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
    SSN = HD_GetSSN; disp(SSN);
    S = LoadSpikesJeff;
    for iCell = 1:length(S.t)
        brainstem_megaplot(iCell)
        set(gcf, 'Position', get(0, 'Screensize'));
        pushdir('D:\Jeff\U01\presentations\GRC conference 2023\megaplots');
        WriteFig(strcat(SSN, '-megaplot', num2str(iCell)))
        disp('figure saved')
        popdir;
    end
    popdir;
end
