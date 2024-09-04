function [nasal_MOVING_amp, nasal_REST_amp, temporal_MOVING_amp, temporal_REST_amp] = calcSaccadeAmplitudeSingleSession(varargin)
%calcSaccadeAmplitude.m
% 2022-10-10. JJS.
% This function tabulates the amplitude of MOVING (evoked) saccades and calculates the amplitude of STATIONARY (spontaneous) saccades for a single session.
%   Inputs
%   Outputs

pixel_ceiling = 45; % *** need to change this 
doPlot = 0;
process_varargin(varargin);

SSN = HD_GetSSN; disp(SSN);

if exist(strcat(SSN, '-saccades-edited.mat')) == 2
    load(strcat(SSN, '-saccades-edited.mat'));
    disp(num2str(k))
    
    [nasal_indices_REST, temporal_indices_REST, ~, ~] = isolateStationarySaccades();
    [~, ~, nasal_indices_MOVING, temporal_indices_MOVING, ~, ~] = isolateManualSaccades();
    
    %     assert(length(nasal_indices_REST) + length(nasal_indices_MOVING) == length(nasalAmplitudes));
    %     assert(length(temporal_indices_REST) + length(temporal_indices_MOVING) == length(temporalAmplitudes));
    
    if length(nasal_indices_REST) + length(nasal_indices_MOVING) ~= length(nasalAmplitudes)
        warning('discrepancy in number of nasal saccades')
    end
    
    if length(temporal_indices_REST) + length(temporal_indices_MOVING) ~= length(temporalAmplitudes)
        warning('discrepancy in number of temporal saccades')
    end
    nasalAmplitudes = nasalAmplitudes(~isnan(nasalAmplitudes)); nasalAmplitudes = nasalAmplitudes*k;
    temporalAmplitudes = temporalAmplitudes(~isnan(temporalAmplitudes)); temporalAmplitudes = temporalAmplitudes*k;
    % JJS. 2024-09-04. Remember that we are now using a conversion factor k for converting from pixels to degrees. 
    nasal_MOVING_amp = nasalAmplitudes(nasal_indices_MOVING);
    nasal_MOVING_amp = nasal_MOVING_amp(abs(nasal_MOVING_amp) < pixel_ceiling);
    temporal_MOVING_amp = temporalAmplitudes(temporal_indices_MOVING);
    temporal_MOVING_amp = temporal_MOVING_amp(abs(temporal_MOVING_amp) < pixel_ceiling);
    
    nasal_REST_amp = nasalAmplitudes(nasal_indices_REST);
    nasal_REST_amp = nasal_REST_amp(abs(nasal_REST_amp) < pixel_ceiling);
    temporal_REST_amp = temporalAmplitudes(temporal_indices_REST);
    temporal_REST_amp = temporal_REST_amp(abs(temporal_REST_amp) < pixel_ceiling);
else
    disp('no EYE data for this session')
    nasal_MOVING_amp = [];
    temporal_MOVING_amp = [];
    nasal_REST_amp = [];
    temporal_REST_amp = [];
end

if doPlot == 1
    bar([median(nasal_MOVING_amp) median(nasal_REST_amp) median(temporal_MOVING_amp) median(temporal_REST_amp)]); hold on
    set(gca, 'XTickLabel', {'nasal M'; 'nasal S'; 'temporal M'; 'temporal S'}, 'FontSize', 20)
    ylabel('saccade amplitude', 'FontSize', 20)
    title(SSN)
    errorbar(1, nanmean(nasal_MOVING_amp), nanstderr(nasal_MOVING_amp), 'Color', 'k')
    errorbar(2, nanmean(nasal_REST_amp), nanstderr(nasal_REST_amp), 'Color', 'k')
    errorbar(3, nanmean(temporal_MOVING_amp), nanstderr(nasal_MOVING_amp), 'Color', 'k')
    errorbar(4, nanmean(temporal_REST_amp), nanstderr(nasal_REST_amp), 'Color', 'k')
end
