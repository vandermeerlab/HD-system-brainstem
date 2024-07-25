function [k, tsdHdeg, diffHdeg] = calibrate_eye_position_single_session_ver2(sd, cfg_in)
% JJS. 2024-05-28.
%   This function attempts to caclulate and plot the HD and eye position traces over one another in order to find the calibration factor for converting
% between pixels and degrees of visual angle. This calculation is based on instructions from Kathy Cullen (see below).
% This version iteratively plots the two traces and asks the user to enter in a coefficient value until he/she is satisfied with the overlap btwn the two traces.
% The conversion to degrees uses the fact that AHV and eye velocity change in step with each other (but in opposite directions). Kathy C. says that w3e
% can approximate degrees by finding the scaling factor at which 0.9*AHV approximate equals the eye velocity trace. This scaling factor usually lies around 0.4 to 06,
% but can vary a lot depending on the particular session. Therefore, it is necessary to estimate k for every single session.
% Future protocols for the task may employ some kind of calibration step, in order to empirically relate the eye position too degrees.
% ***Note to self: ask Brandie Verdone how she does the calibration step in the Cullen lab.

% JJS. 2024-07-25. ver2
%   This version additionally saves the factor k and the new eye variables in deg/sec into the *saccades-edited.mat file.
%   The names for the new variables are tsdHdeg.data; diffHdeg.data

% 1) differentiate the eye movement trace in pixels (accounting for time) and then
% 2) compare this to the head velocity that was applied - hopefully you recorded head position (fingers crossed)? You should convert this head poistion signal to degrees and then differentiation this (again accounting for time).
% 3) then focus on the slow phase segments in between the quick phases
% 4) then find a calibration factor such that the eye velocity trace is on average 0.9 of head velocity
%   Inputs
%
%   Outputs
%           k - conversion factor for the session. This is the coefficient for eye velocity in order to get to ~0.9 head velocity (see Kathy's email above)
%           tsdHdeg - 1x1 struct with fields
%                   .type - 'tsd'
%                   .units - 'degrees' [this is pupil position in the eye orbit]
%                   .tvec - timebase
%                   .data - data
%                   .label - cell with string, indicating the meaning of this variable
%           diffHdeg - 1x1 struct with fields same as above ^^
%                   .unit - 'degrees per second' [this is pupil velocity]

cfg_def.doSave = 1;
cfg_def.FontSize = 25;
cfg_def.MarkerSize = 15;
cfg_def.doSmooth = 1;
cfg_out = ProcessConfig(cfg_def,cfg_in);

% Initialize
SSN = HD_GetSSN;
k = [];
tsdHdeg = [];  % this is a tsd of the eye position in degrees (as approximated by the 0.9 AHV rule)
diffHdeg = []; % this is a tsd of eye velocity in degrees per second

%% Make sure that eye data is present
filename = strcat(SSN, '-saccades-edited.mat');
doExist = exist(filename);
if doExist ~= 2
    warning('-saccades-edited.mat file could not be found!')
    return
else
    load(filename)
    disp('saccade data loaded')
end

%% Isolate Slow Phase epochs
cfg_slow.doPlot = 0;
disp('restricting to Slow Phase only')
[data_out, data_outR, nSpikesRemoved] = getSlowPhaseData(cfg_slow, sd, []); % iCell is what neuron to look at. We don't care about the neural data here, but it is a mandatory field.
%  data_outR.____                               dataoutR - a structure with the 'de-saccaded' (restricted) data
%          horiz_eye_pos: [1×1 struct]
%          horiz_eye_vel: [1×1 struct]
%   horiz_eye_vel_smooth: [1×1 struct]
%                    AHV: [1×1 struct]
%                      S: [1×1 struct]
%                     HD: [1×1 struct]
%                     fr: [1×1 struct]
A = abs(data_outR.horiz_eye_vel.data) > 200; % remove abnormally high or low (artifactual) values for the eye velocity trace
data_outR.horiz_eye_vel.data(A) = NaN; % make those values NaN
data_outR.AHV.data(A) = NaN; % make the corresponding AHV values NaN

if cfg_out.doSmooth
    eyeR = tsd(data_outR.horiz_eye_vel.tvec, smoothdata(data_outR.horiz_eye_vel.data));
    ahvR = tsd(data_outR.AHV.tvec, smoothdata(-data_outR.AHV.data)); % introducing a sign change here so that AHV and EYE are correlated.
    disp('using smoothed data')
else
    eyeR = tsd(data_outR.horiz_eye_vel.tvec, data_outR.horiz_eye_vel.data);
    ahvR = tsd(data_outR.AHV.tvec, -data_outR.AHV.data); % introducing a sign change here so that AHV and EYE are correlated.
end

k = 0.40;
cont = 1; % 'continue' variable. Keep running the for loop until user is satisfied with the coefficient
while cont == 1
    clf; hold on
    plot(ahvR.tvec, (0.9).*ahvR.data, '-o', 'Color', 'r', 'MarkerSize', cfg_out.MarkerSize)
    plot(eyeR.tvec, k.*eyeR.data, '-o', 'Color', 'b', 'MarkerSize', cfg_out.MarkerSize)
    %     plot(eye.tvec, k.*eye.data, 'b')
    set(gca, 'FontSize', cfg_out.FontSize)
    % c = axis;
    %     axis([177 192 -250 250]);
    if cfg_out.doSmooth
    title(strcat(sd.SSN, ' smoothed'))
    else
        title(sd.SSN)
    end
    legend('9/10*AHV (deg/s)', strcat(num2str(k),'*EV (pixels/s)'), 'Location', 'NorthWest')
    x = input('Are you satisfied with the overlap? Type y or n    ', "s");
    if strcmp(x, 'n') == 1
        y = input('Type in a new numerical value for the coefficient     ');
        assert(isnumeric(y))
        k = y;
    elseif strcmp(x, 'y') == 1
        disp('k value accepted')
        cont = 0;
    else
        error('must type y or n as a reply')
    end
end

tsdHdeg = tsd(tsdH.tvec, tsdH.data.*k); % tsdH is the pupil position, in units of pixels
tsdHdeg.units = 'degrees';
tsdHdeg.label = {'calibrated pupil position'};
diffHdeg = tsd(diffH.tvec, diffH.data.*k); % diffH is the pupil velocity, in units of pixels per second
diffHdeg.units = 'degrees per second';
diffHdeg.label = {'calibrated pupil velocity'};

if cfg_out.doSave
    save(filename, 'k', 'tsdHdeg', 'diffHdeg', "-append")
    disp('data saved to saccades-edited.mat file')
    savefig(strcat(SSN, '_pupil velocity calibrated'))
end


