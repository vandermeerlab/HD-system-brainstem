
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
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}))
    SSN = HD_GetSSN; disp(SSN);
    [temporalSaccades(iSess), nasalSaccades(iSess), combinedSaccades(iSess), index_tP_final, index_nP_final, tsdH, tsdV, diffH, diffV] = processPupilData2([]);
    popdir;
end


