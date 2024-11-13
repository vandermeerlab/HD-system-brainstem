function [wheel_speed] = get_wheel_speed_tsd_only(varargin)

process_varargin(varargin); 

sd = LoadSessionData([]);
updownTSD = getQEupdown([]);
state_tsd = ConvertQEUpDownToState(updownTSD);
[~, wheel_tsd] = ConvertQEStatesToAngle([], state_tsd);
[d, speed, cfg_w] = ConvertWheeltoSpeed([], wheel_tsd);   % d = distance.  speed = speed of the wheel in cm/sec
wheel_speed.tvec = speed.tvec - sd.starttime;
wheel_speed.data = -speed.data; % we want forward motion to be displayed as a positive velocity


