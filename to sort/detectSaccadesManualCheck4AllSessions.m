function detectSaccadesManualCheck4AllSessions(cfg_in)
cd('C:\Jeff\U01\datatouse');
fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    detectSaccadesManualCheck4(cfg_in);
    popdir;
end