function [ID, LRatio, S_N_R] = collect_all_cluQual()
format long
fd = FindFiles('*keys.m');
ft = FindFiles('*.t');
ID = nan(1,length(ft));
LRatio = nan(1,length(ft));
S_N_R = nan(length(ft), 4);
iCell = 0;
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN; disp(SSN);
    [fc] = FindFiles(strcat(SSN, '*.t'));
    [~, b, ~] = fileparts(fc);
    if iscell(b)
        cellID = b;
    else
        cellID = {b};
    end
    otherID = cellID;
    for iZ = 1:length(cellID)
        newStr = strrep(otherID{iZ} ,'_', '-');
        %         otherID{iZ}(21) = '-';
        iCell = iCell + 1;
        if exist(strcat(newStr, '-ClusterQual.mat'))
            load(strcat(newStr, '-ClusterQual.mat'));
            ID(iCell) = CluSep.IsolationDist;
            LRatio(iCell) = CluSep.Lratio;
            S_N_R(iCell,:) = CluSep.SNR;
        else
            disp('no clusterqual file found!')
        end
    end
end