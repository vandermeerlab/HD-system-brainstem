function [labels, times, fdr] = getEventLabels_Allsessions(fd, startSess, endSess)
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
    [a, b, c] = fileparts(fd{iSess});
    fdr{iSess} = b(1:end-5);
    pushdir(fileparts(fd{iSess}));
    events_ts = LoadEvents([]);
    SSN = HD_GetSSN; disp(SSN);
    for iLabel = 1: length(events_ts.label)
        labels{iSess,iLabel} = events_ts.label{iLabel};
        times{iSess,iLabel} = events_ts.t{iLabel};
    end
end
fdr = fdr';