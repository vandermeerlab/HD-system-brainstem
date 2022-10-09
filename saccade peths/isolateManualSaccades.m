function [Mstart, Mend, nasal_indices_MOVING, temporal_indices_MOVING, nasal_timestamps_MOVING, temporal_timestamps_MOVING] = isolateManualSaccades
%JJS. 2022-09-12. Identify start and end times for when the platform is moving.
%   Get timestamps for (evoked) saccades that occur when the platform is moving.
%Mstart = start times for platform-moving (M) periods;
%Mend = end times for platform-moving (M) periods.


% JJS. 2022-09-01.
% function to pull out only saccades from periods where the recording platform was stationary

% inputs: start and end times for stationary epoch [from *-AHV_StationaryTimes.mat]
% outputs:  nasal_indices_MOVING:           nasal saccade indices (mouse stationary)
%           temporal indices:               temporal saccade indices (mouse stationary)
%           nasalREST:                      timestamps for nasal saccades (mouse stationary)
%           temporalREST:                   timestamps for temporal saccades (mouse stationary)
nasal_indices_MOVING = [];
temporal_indices_MOVING = [];


SSN = HD_GetSSN; disp(SSN);
if exist(strcat(SSN, '-AHV_StationaryTimes.mat'))==2 && exist(strcat(SSN, '-saccades-edited.mat')) ==2
    load(strcat(SSN, '-AHV_StationaryTimes.mat'));
    load(strcat(SSN, '-saccades-edited.mat'));
    assert(length(STtstart) == length(STtend))
    %     nasal_indices_MOVING = [];
    %     temporal_indices_MOVING = [];
    
    nasalSaccadesToUse = nasalSaccades(~isnan(nasalSaccades));
    temporalSaccadesToUse = temporalSaccades(~isnan(temporalSaccades));
    
%     events_ts = LoadEvents([]);
    
    Mstart = STtend(1:end-1);
    Mend = STtstart(2:end);
    
    assert(length(Mstart) == length(Mend))
    
    ST = length(Mstart); % number of stationary periods
    for iST = 1:ST
        overlap_N = nasalSaccadesToUse > Mstart(iST) & nasalSaccadesToUse < Mend(iST);
        N = find(overlap_N);
        nasal_indices_MOVING = horzcat(nasal_indices_MOVING, N);
    
        overlap_T = temporalSaccadesToUse > Mstart(iST) & temporalSaccadesToUse < Mend(iST);
        T = find(overlap_T);
        temporal_indices_MOVING = horzcat(temporal_indices_MOVING, T);
        disp(iST)
    end
    
    nasal_timestamps_MOVING = nasalSaccadesToUse(nasal_indices_MOVING);
    temporal_timestamps_MOVING = temporalSaccadesToUse(temporal_indices_MOVING);
    
    
    
else
    warning('one or more files missing')
    Mstart = NaN;
    Mend = NaN;
end

