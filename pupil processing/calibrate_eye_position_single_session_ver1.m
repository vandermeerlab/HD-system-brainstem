function [C] = calibrate_eye_position_single_session_ver1(sd, cfg_in) 
% JJS. 2024-05-28. 
% This function attempts to caclulate and plot the HD and eye position traces over one another in order to find the calibration factor for converting 
% between pixels and degrees of visual angle. This calculation is based on instructions from Kathy Cullen (see below). 

% 1) differentiate the eye movement trace in pixels (accounting for time) and then 
% 2) compare this to the head velocity that was applied - hopefully you recorded head position (fingers crossed)? You should convert this head poistion signal to degrees and then differentiation this (again accounting for time).
% 3) then focus on the slow phase segments in between the quick phases 
% 4) then find a calibration factor such that the eye velocity trace is on average 0.9 of head velocity
%   Inputs  
%                   
%   Outputs
%           C - conversion factor for the session. This is the coefficient for eye velocity in order to get to ~0.9 head velocity (see Kathy's email above) 

cfg_def = [];
cfg_out = ProcessConfig(cfg_def,cfg_in);
%% Isolate Slow Phase epochs
cfg_slow.doPlot = 0;
[data_out, ~] = getSlowPhaseData(cfg_slow, sd, 1); % iCell is what neuron to look at. We don't care about the neural data here, but it is a mandatory field.
 %          horiz_eye_pos: [1×1 struct]
 %          horiz_eye_vel: [1×1 struct]
 %   horiz_eye_vel_smooth: [1×1 struct]
 %                    AHV: [1×1 struct]
 %                      S: [1×1 struct]
 %                     HD: [1×1 struct]
 %                     fr: [1×1 struct]




























%           sd - session data structure for session of interest. Make sure you are in the session folder. Use sd = LoadSessionData([]);
%                   sd.AHV - angular head velocity (in units of degrees / second)
%                   sd.diffH - horizontal eye velocity (in pixels / second)
