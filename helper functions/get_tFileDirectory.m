function [tt_name] = get_tFileDirectory(fd)

tt_name = {};
cellcounter = 0;
if isempty(fd)
    fd = FindFiles('*keys.m');
end
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}))
    SSN = HD_GetSSN; disp(SSN);
    spikefiles = FindFiles('*.t');
    [a b c] = fileparts(spikefiles);
    %     ind = strfind(b, 'TT');
    for iCell = 1:length(spikefiles);
        cellcounter = cellcounter + 1;
        if iscell(b)
            tt_name{cellcounter,1} = b{iCell};
        else
            tt_name{cellcounter,1} = b;
        end
        disp(num2str(cellcounter))
    end
    popdir;
end








