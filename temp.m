% for iSess = 1:length(fd)
%     [~, b, ~] = fileparts(fd{iSess});
%     fdr{iSess} = b;
% end

for iCell = 1:length(Z)
    percent_n(iCell) = Z{iCell}.percent_n(1,end);
    percent_t(iCell) = Z{iCell}.percent_t(1,end);
end