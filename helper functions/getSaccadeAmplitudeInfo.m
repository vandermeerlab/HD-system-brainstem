function [temporal_mean, nasal_mean, diff_amp] = getSaccadeAmplitudeInfo(fd, cfg_in)
% JJS. 2024-07-29.
% This function gathers the average temporal and nasal saccade amplitdues for each session, as well as the difference between them. This info. is 
% of interest because we want to look to see if there are any systematic differences in saccade amplitude. 
% ***Note: these values are in pixels. Should rewrite to use degrees of visual angle.

% Inputs: 
%           fd - file directory list. if empty, will use current folder
%           cfg - config variables
% Outputs
%           temporal mean - mean of temporal saccade amplitudes
%           nasal mean - mean of nasal saccades amplitdues
%           diff_amp - difference in (absolute) amplitudes of each saccade type 
                    % A NEGATIVE value for diff_amp means that NASAL saccades were larger. 
                    % A POSITIVE value for diff_amp means that TEMPORAL saccades were larger. 

if isempty(fd)
    fd = FindFiles('*keys.m');
end

cfg_def.startSess = 1;
cfg_def.endSess = length(fd);
cfg_out = ProcessConfig(cfg_def,cfg_in);

for iSess = cfg_out.startSess : cfg_out.endSess
    pushdir(fileparts(fd{iSess}));
    %% Make sure that eye tracking exists for this session
    SSN = HD_GetSSN; disp(SSN)
    if exist(strcat(SSN, '-saccades-edited.mat')) == 2 
        load(strcat(SSN, '-saccades-edited.mat'))
        temporal_mean(iSess) = nanmean(temporalAmplitudes); % mean temporal saccade amplitude for this session
        nasal_mean(iSess) = nanmean(nasalAmplitudes); % mean nassal saccade amplitdue for this session 
        diff_amp(iSess) = abs(temporal_mean(iSess)) - abs(nasal_mean(iSess)); % difference of the average amplitudes for this session 
    end
    popdir;
    clear temporalAmplitudes
    clear nasalAmplitudes 
end