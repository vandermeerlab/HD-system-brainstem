for iSess = 1:length(fd)
    [~, b, ~] = fileparts(fd{iSess});
    fdr{iSess} = b;
end