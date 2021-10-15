function [Z] = saccade_peth_shuffled_saccade_times_all_sessions(cfg_in)
% 10/2021. JJS.
% get Zscore value for nasal and temporal saccades for a given cell. Uses some -pre and -post saccade window for consideration. Loops over sessions. 

cfg_def = [];
cfg = ProcessConfig(cfg_def,cfg_in); 

fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));     
    SSN = HD_GetSSN; disp(SSN); 

    cfg_pupil = [];
    cfg_pupil.threshAdj  = 4;
    cfg_pupil.threshH = 10;
    cfg_pupil.threshL = -10;
    cfg_pupil.scalingfactor = 1;  % for shrinking the pupil trace so its the same height as diffH
    cfg_pupil.artifactThresh = 4;  % units of pixels sq.
    [temporalSaccades, nasalSaccades, combinedSaccades, index_tP_final, index_nP_final, tsdH, tsdV, diffH, diffV] = processPupilData2(cfg_pupil, 'doPlotThresholds', 0, 'doPlotEverything', 0);
    
    S = LoadSpikes([]);
    [M, Mfr] = saccade_peth_shuffled_saccade_times(myCell, t); % Creates a matrix of shuffled spike times. Takes as input a spike train (ts) and saccade times (t)
    
    
end

    
    
    