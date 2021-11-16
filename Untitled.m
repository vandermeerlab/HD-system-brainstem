title('This_title has an underline', 'Interpreter', 'none');



cellcounter = 0;
fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}))
    SSN = HD_GetSSN; disp(SSN);
    spikefiles = FindFiles('*.t');
    [a b c] = fileparts(spikefiles);
    ind = strfind(b, 'TT');
    for iCell = 1:length(ind);
        cellcounter = cellcounter + 1;
        tt_place(iCell) = b{iCell}(ind{iCell}+3);
    end
    popdir;
end
arrayfun(@(s) set(s,'EdgeColor','none'), findobj(gcf,'type','surface'))


fd = FindFiles('*keys.m');
temporalSaccades = nan(1,length(fd));
nasalSaccades = nan(1,length(fd));
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}))
    SSN = HD_GetSSN; disp(SSN);
    sd = LoadSessionData([], 'EYE', false);
    if exist(strcat(SSN, '-VT1.nvt')) == 2;
    [temporalSaccades(iSess), nasalSaccades(iSess), combinedSaccades(iSess), index_tP_final, index_nP_final, tsdH, tsdV, diffH, diffV] = processPupilData2([]);
    else
        disp('skipping')
    end
    popdir;
end


