% counts up how many nasal and temporal saccades are present in each session
cd('C:\Jeff\U01\datatouse');
fd = FindFiles('*keys.m');
num_temporalSaccades = nan(1,length(fd));
num_nasalSaccades = nan(1,length(fd));
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}))
    SSN = HD_GetSSN; disp(SSN);
    sd = LoadSessionData([], 'EYE', false);
    if exist(strcat(SSN, '-VT1.smi')) == 2
        cfg_in = [];
        cfg_in.doPlotThresholds = 1;
        cfg_in.doPlotEverything = 1;
        [temporalSaccades, nasalSaccades, ~, ~, ~, ~, ~, ~, ~] = processPupilData2(cfg_in, sd);
        num_temporalSaccades(iSess) = length(temporalSaccades);
        num_nasalSaccades(iSess) = length(nasalSaccades);
        
    else
        disp('skipping')
    end
    popdir;
end