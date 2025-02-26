function classify_saccades_PFD(cfg_in, sd)
% 2025-01-28. JJS. This function labels saccades as TOWARD or AWAY from the preferred firing direciton (PFD) of the neuron, if it is a head direction cell.
%   There are several options for being more (or less) restrictive with the data, including AHV speed cutoff, saccade amplitude threshold, how close the mouse needs
%   to be to the PFD in order to use the saccade, definition of IN-FIELD saccades vs. OUT-OF-FIELD saccades, and splitting things up between VOR saccades and
%   "voluntary" saccades (when the mouse is stationary, reacting to a sound or light stimulus, or simply choosing to saccade).

% Inputs:
%           sd          - session data. A structure that includes ...  sd.AHV
%           cfg_in      - input parameters

% Outputs:
%           TC
%


cfg_def.AHV_min = 1;
cfg_def.AHVthresh = 20; % cm/sec

cfg = ProcessConfig(cfg_def,cfg_in);
SSN = HD_GetSSN;
load(strcat(SSN, '-saccades-edited.mat'), 'nasalSaccades', 'temporalSaccades', 'combinedSaccades')

%% Restrict to low AHV times, if desired
if cfg.AHV_min == 1
    % How much of session is ROTATION
    total_AHV_samples = length(sd.AHV.tvec);
    rotation_time = antirestrict(sd.AHV, sd.STstart, sd.STend);
    rotation_AHV_samples = length(rotation_time.tvec);
    fraction_rotation = rotation_AHV_samples/total_AHV_samples;
    disp(strcat('fraction rotation =', num2str(round(fraction_rotation,2))))
    % How much of session is low AHV?
    low_AHV = abs(sd.AHV.data) < cfg.AHVthresh;
    low_AHV_samples = sum(low_AHV);
    low_AHV_seconds = low_AHV_samples/(1/median(diff(sd.AHV.tvec)));
    low_AHV_minutes = low_AHV_seconds/60;
    fraction_low_AHV = low_AHV_samples / length(sd.AHV.tvec); 
    disp(strcat('minutes of low AHV time =', num2str(round(low_AHV_minutes,2))))
    disp(strcat('fraction low AHV samples =', num2str(round(fraction_low_AHV,2))))
    
    % Find start/stop times for AHV restriction
    low_AHV_diff = horzcat(NaN, diff(low_AHV)); % take diff of the indices to get transition points
    % High to Low
    % ---------------------------------------
    %                                       -
    %                                       ------------------------------------------------
    low_AHV_ones = find(low_AHV_diff == 1); sum_low_AHV_ones = length(low_AHV_ones); % find the transition points from high(er) AHV to low AHV
    low_AHV_tstart = sd.AHV.tvec(low_AHV_ones); 
    low_AHV_tstart = low_AHV_tstart';
    % Low to High
    %                                       ------------------------------------------------
    %                                       -
    % ---------------------------------------
    low_AHV_minus_ones = find(low_AHV_diff == -1); sum_low_AHV_minus_ones = length(low_AHV_minus_ones);  % find transitions from low AHV back to high AHV
    low_AHV_tend = sd.AHV.tvec(low_AHV_minus_ones); 
    low_AHV_tend = low_AHV_tend';
    % Account for the first value being low or high. *** Need to account for all possiblities ***
    if low_AHV_tstart(1) > low_AHV_tend(1)   % In other words, if the first transition is from high to low, then the session started out LOW. Make tstart(1) = the first timestamp
        low_AHV_tstart = [sd.AHV.tvec(1) low_AHV_tstart];
    end
    if low_AHV_tstart(end) > low_AHV_tend(end) % In other words, if the last transition is from high to low, then the session ended LOW. Make tend(end) = last timestamp
        low_AHV_tend = [low_AHV_tend sd.AHV.tvec(end)];
    end
    assert(length(low_AHV_tstart) == length(low_AHV_tend))
    
    
    AHV_tsd = tsd(sd.AHV.tvec(low_AHV), sd.AHV.data(low_AHV));
    
    % How many saccades were removed?
    
    
    
    % Plot full data vs. restricted data
    clf
    plot(sd.AHV.tvec, sd.AHV.data); hold on
    plot(AHV_tsd.tvec, AHV_tsd.data); xlabel('time'); ylabel('AHV')
    disp('AHV (all) and AHV (restricted) are plotted. Press any key to continue')
    
    
    
    
else
    AHV_tsd = sd.AHV;
end
%%

























