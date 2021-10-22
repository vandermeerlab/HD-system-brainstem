function [zT, zN] = saccade_peth_shuffled_saccade_times_all_sessions(cfg_in)
% 10/2021. JJS.
% get Zscore value for nasal and temporal saccades for a given cell. Uses some -pre and -post saccade window for consideration. Loops over sessions.

cfg_def = [];
cfg_def.doPlotShuffles = 0;
cfg = ProcessConfig(cfg_def,cfg_in);

fd = FindFiles('*keys.m');
disp(length(fd))
cellCounter = 0;
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN; disp(SSN);
    if exist(strcat(SSN, '-VT1.smi')) == 2
        
        %% Calculate the saccade times (for now, on the fly)
        cfg_pupil = [];
        cfg_pupil.threshAdj  = 4;
        cfg_pupil.threshH = 10;
        cfg_pupil.threshL = -10;
        cfg_pupil.scalingfactor = 1;  % for shrinking the pupil trace so its the same height as diffH
        cfg_pupil.artifactThresh = 4;  % units of pixels sq.
        disp('calculating saccade times')
        [temporalSaccades, nasalSaccades, ~, ~, ~, ~, ~, ~, ~] = processPupilData2(cfg_pupil, 'doPlotThresholds', 0, 'doPlotEverything', 0);
        
        S = LoadSpikesJeff();
        for iCell = 1:length(S.t)
            cellCounter = cellCounter + 1;
            myCell = SelectTS([], S, iCell);
            [~, ~, zT(cellCounter,:), zN(cellCounter,:)] = saccade_peth_shuffled_saccade_times2(cfg_in, myCell, iCell, temporalSaccades, nasalSaccades);
        end       
    else
        disp('No video timestamps exist for this session. Skipping...')
    end
    popdir;
end


