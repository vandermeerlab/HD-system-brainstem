function [nasal_indices_REST, temporal_indices_REST, nasal_timestamps_REST, temporal_timestamps_REST] = isolateStationarySaccades() 
% JJS. 2022-09-01.
% function to pull out only saccades from periods where the recording platform was stationary

% inputs: start and end times for stationary epoch [from *-AHV_StationaryTimes.mat]
% outputs:  STtstart:                           start of stationary periods (timestamps)
%           STtend:                             end of stationary periods (timestamps)
%           nasal_indices_REST:                 nasal saccade indices (mouse stationary)
%           temporal indices:                   temporal saccade indices (mouse stationary)
%           nasal_timestamps_REST:              timestamps for nasal saccades (mouse stationary)
%           temporal_timestamps_REST:           timestamps for temporal saccades (mouse stationary)

SSN = HD_GetSSN; % disp(SSN);
if exist(strcat(SSN, '-AHV_StationaryTimes.mat'))==2 && exist(strcat(SSN, '-saccades-edited.mat')) ==2
    load(strcat(SSN, '-AHV_StationaryTimes.mat'));
    load(strcat(SSN, '-saccades-edited.mat'));
    assert(length(STtstart) == length(STtend))
    nasal_indices_REST = [];
    temporal_indices_REST = [];
    
    nasalSaccadesToUse = nasalSaccades(~isnan(nasalSaccades));
    temporalSaccadesToUse = temporalSaccades(~isnan(temporalSaccades));
    
    ST = length(STtstart); % number of stationary periods
    for iST = 1:ST
        overlap_N = nasalSaccadesToUse > STtstart(iST) & nasalSaccadesToUse < STtend(iST);
        N = find(overlap_N);
        nasal_indices_REST = horzcat(nasal_indices_REST, N);
        
        overlap_T = temporalSaccadesToUse > STtstart(iST) & temporalSaccadesToUse < STtend(iST);
        T = find(overlap_T);
        temporal_indices_REST = horzcat(temporal_indices_REST, T);
        % disp(iST)
    end
    
    nasal_timestamps_REST = nasalSaccadesToUse(nasal_indices_REST);
    temporal_timestamps_REST = temporalSaccadesToUse(temporal_indices_REST);
    
else
    warning('one or more files missing')
    nasal_indices_REST = NaN;
    temporal_indices_REST = NaN;
    nasal_timestamps_REST = NaN;
    temporal_timestamps_REST = NaN;
end



