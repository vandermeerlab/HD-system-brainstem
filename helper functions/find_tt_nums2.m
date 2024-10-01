function [tt_place] = find_tt_nums2(tfilelist)
% JJS. 2024-10-01. This function searches through the neurons in tfilelist to determine how many cells were recorded from each tetrode (over all cells)

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








