function [wheel] = get_wheel_speed

% cfg = ProcessConfig(cfg_def, cfg_in);

updownTSD = getQEupdown([]);
state_tsd = ConvertQEUpDownToState(updownTSD);
[~, wheel_tsd] = ConvertQEStatesToAngle([], state_tsd);
[d, speed, cfg_w] = ConvertWheeltoSpeed([], wheel_tsd);
wheel_speed = -speed.data; % we want forward motion to be displayed as a positive velocity
wheel = tsd(speed.tvec, wheel_speed);