function [wheel] = get_wheel_speed(starttimeUNIX)

%   Inputs:
%               starttimeUNIX - timestamp at which the session started, in UNIX time. This can be obtained from the encoder CSC or an LFP CSC. 
%   Outputs:
%               wheel - TSD of wheel speed in centimeters per second, starting at t = 0 seconds. 
% cfg = ProcessConfig(cfg_def, cfg_in);

updownTSD = getQEupdown([]);
state_tsd = ConvertQEUpDownToState(updownTSD);
[~, wheel_tsd] = ConvertQEStatesToAngle([], state_tsd);
[d, speed, cfg_w] = ConvertWheeltoSpeed([], wheel_tsd);
wheel_speed = -speed.data; % we want forward motion to be displayed as a positive velocity
speed.tvec = speed.tvec - starttimeUNIX; 
wheel = tsd(speed.tvec, wheel_speed);
