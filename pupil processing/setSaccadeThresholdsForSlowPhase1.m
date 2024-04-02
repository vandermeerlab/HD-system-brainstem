function setSaccadeThresholdsForSlowPhase1(cfg_in)
% JJS. 2024-03-28. Function for manually scrolling through each session in the directory and setting a threshold (in pixels) for temporal and nasal saccades. 
% A set of temporal and nasal saccade times are then calculated and saved to a .mat file in the session folder. These thresholds are specifically for slow phase 
% eye movement analyses (tuning curves). For quick phase analyses / analyses of the unrestricted data, use the "SSN-YYYY-MM-DD-saccades-edited.mat file."

fd = FindFiles('*keys.m');

for iSess = 1:length(fd)
    disp(iSess)
    pushdir(fileparts(fd{iSess}));
    sd = LoadSessionData([]);
    
    
    
    
    
    
    popdir;
end

   