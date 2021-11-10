% counts up how many nasal and temporal saccades are present in each session
cd('C:\Jeff\U01\datatouse');
fd = FindFiles('*keys.m');
num_temporalSaccades = nan(1,length(fd));
num_nasalSaccades = nan(1,length(fd));
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}))
    SSN = HD_GetSSN; disp(SSN);
    if exist(strcat(SSN, '-VT1.smi')) == 2
        cfg = [];
        cfg.doPlotThresholds = 1;
        cfg.doPlotEverything = 1;
        [temporalSaccades, nasalSaccades, ~, ~, ~, ~, ~, ~, ~] = processPupilData2(cfg);
        num_temporalSaccades = length(temporalSaccades);
        num_nasalSaccades = length(nasalSaccades);
        
    else
        disp('skipping')
    end
    popdir;
end