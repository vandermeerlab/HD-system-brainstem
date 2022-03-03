updownTSD = getQEupdown([]);
state_tsd = ConvertQEUpDownToState(updownTSD);
[angle_tsd, wheel_tsd] = ConvertQEStatesToAngle([], state_tsd);
[d, speed, cfg] = ConvertWheeltoSpeed([], wheel_tsd);