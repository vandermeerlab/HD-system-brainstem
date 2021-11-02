
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
        if iscell(b)
            tt_place(cellcounter) = str2num(b{iCell}(ind{iCell}+3));
        else
            tt_place(cellcounter) = str2num(b(ind+3));
        end
    end
    popdir;
end








