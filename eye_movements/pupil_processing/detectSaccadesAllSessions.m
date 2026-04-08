fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN; disp(SSN);
    if exist(strcat(SSN, '-VT1_proc.mat')) ==2 
        disp('video file found')
        [m] = detectSaccadesManualCheck([]);
    else
        disp('video file not detected. Skipping.')
    end
    popdir;
end
