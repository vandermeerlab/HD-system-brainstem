function [k] = calibrate_eye_position_single_session_ver1(sd, cfg_in)
% JJS. 2024-05-28.
% This function attempts to caclulate and plot the HD and eye position traces over one another in order to find the calibration factor for converting
% between pixels and degrees of visual angle. This calculation is based on instructions from Kathy Cullen (see below).
% This version iteratively plots the two traces and asks the user to enter in a coefficient value until he/she is satisfied with the overlap btwn the two traces.

% 1) differentiate the eye movement trace in pixels (accounting for time) and then
% 2) compare this to the head velocity that was applied - hopefully you recorded head position (fingers crossed)? You should convert this head poistion signal to degrees and then differentiation this (again accounting for time).
% 3) then focus on the slow phase segments in between the quick phases
% 4) then find a calibration factor such that the eye velocity trace is on average 0.9 of head velocity
%   Inputs
%
%   Outputs
%           C - conversion factor for the session. This is the coefficient for eye velocity in order to get to ~0.9 head velocity (see Kathy's email above)
cfg_def.FontSize = 25;
cfg_def.MarkerSize = 15;
cfg_def.doSmooth = 1;
cfg_out = ProcessConfig(cfg_def,cfg_in);
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
eyeR = tsd(data_outR.horiz_eye_vel.tvec, data_outR.horiz_eye_vel.data);
ahvR = tsd(data_outR.AHV.tvec, -data_outR.AHV.data); % introducing a sign change here so that AHV and EYE are correlated.
k = 0.40;
cont = 1; % 'continue' variable. Keep running the for loop until user is satisfied with the coefficient
while cont == 1
    clf; hold on
    plot(ahvR.tvec, (0.9).*ahvR.data, '-o', 'Color', 'r', 'MarkerSize', cfg_out.MarkerSize)
    plot(eyeR.tvec, k.*eyeR.data, '-o', 'Color', 'b', 'MarkerSize', cfg_out.MarkerSize)
%     plot(eye.tvec, k.*eye.data, 'b.', 'MarkerSize', cfg_out.MarkerSize)
%     plot(ahv.tvec, (0.9).*ahv.data, 'r', 'MarkerSize', cfg_out.MarkerSize)

%     plot(ahv.tvec, (0.9).*ahv.data, 'r.', 'MarkerSize', cfg_out.MarkerSize)
%     plot(eye.tvec, k.*eye.data, 'b.', 'MarkerSize', cfg_out.MarkerSize)
%     plot(eye.tvec, k.*eye.data, 'b.', 'MarkerSize', cfg_out.MarkerSize)
%     plot(ahv.tvec, (0.9).*ahv.data, 'r', 'MarkerSize', cfg_out.MarkerSize)
%     plot(eye.tvec, k.*eye.data, 'b')
    set(gca, 'FontSize', cfg_out.FontSize)
    % c = axis;
    %     axis([177 192 -250 250]);
    title(sd.SSN)
    legend('9/10*AHV (deg/s)', strcat(num2str(k),'*EV (pixels/s)'), 'Location', 'NorthWest')
    x = input('Are you satisfied with the overlap? Type y or n    ', "s");
    if strcmp(x, 'y') == 1
        return
    elseif strcmp(x, 'n') == 1
        y = input('Type in a new numerical value for the coefficient     ');
        assert(isnumeric(y))
        k = y;
    else
        error('must type y or n as a reply')
    end
end
