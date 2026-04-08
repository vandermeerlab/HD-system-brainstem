fd = FindFiles('*keys.m');
for iSess = 1:length(fd);
    [~, fdr{iSess}, ~] = fileparts(fd{iSess});
end
fdr = fdr';