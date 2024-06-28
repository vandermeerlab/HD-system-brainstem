for iSess = 1:length(fd);
    pushdir(fileparts(fd{iSess}))
    [a{iSess}, b{iSess}, c{iSess}] = fileparts(fd{iSess});
end
