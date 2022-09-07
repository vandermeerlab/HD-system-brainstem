function [nasal_indices, temporal_indices, nasalREST, temporalREST] = isolateStationarySaccades()
% JJS. 2022-09-01.
% function to pull out only saccades from periods where the recording platform was stationary

% inputs: start and end times for stationary epoch [from *-AHV_StationaryTimes.mat]
% outputs:  nasal_indices:          nasal saccade indices (mouse stationary)
%           temporal indices:       temporal saccade indices (mouse stationary)
%           nasalREST:              timestamps for nasal saccades (mouse stationary)
%           temporalREST:           timestamps for temporal saccades (mouse stationary) 

SSN = HD_GetSSN; disp(SSN);
if exist(strcat(SSN, '-AHV_StationaryTimes.mat'))==2 && exist(strcat(SSN, '-saccades-edited.mat')) ==2
    load(strcat(SSN, '-AHV_StationaryTimes.mat'));
    load(strcat(SSN, '-saccades-edited.mat')); 
    assert(length(STtstart) == length(STtend))
    nasal_indices = [];
    temporal_indices = [];
    
    nasalSaccadesToUse = nasalSaccades(~isnan(nasalSaccades));
    temporalSaccadesToUse = temporalSaccades(~isnan(temporalSaccades)); 
    
    ST = length(STtstart); % number of stationary periods
    for iST = 1:ST
        overlap_N = nasalSaccadesToUse > STtstart(iST) & nasalSaccadesToUse < STtend(iST);
        N = find(overlap_N);
        nasal_indices = horzcat(nasal_indices, N);
        
        overlap_T = temporalSaccadesToUse > STtstart(iST) & temporalSaccadesToUse < STtend(iST);
        T = find(overlap_T);
        temporal_indices = horzcat(temporal_indices, T);
        disp(iST)
    end
    
    nasalREST = nasalSaccadesToUse(nasal_indices);
    temporalREST = temporalSaccadesToUse(temporal_indices); 

else
    warning('one or more files missing')
    nasal_indices = NaN;
    temporal_indices = NaN;
    nasalREST = NaN; 
    temporalREST = NaN; 
end



