function [angle_tsd, wheel_tsd] = ConvertQEStatesToAngle(cfg_in, state_tsd)


cfg_def = [];
cfg_def.dphi = 1; % encoder angle step -- property of the encoder (should be obtained by calibration procedure)
cfg_def.angle_to_turn  = 1024; % This is how many pulses per revolution for the wheel encoder. We have a Grayhill encoder from digikey. The part number is 63KS256. [this (angle_to_turn) is the conversion factor for coverting integer "angle" values into # of wheel turns. This was determined empirically with a recording session on 3/3/2022.] 


% state definitions:
% high-high: 1
% high-low: 2
% low-high: 3
% low-low: 4

% this variable maps state transitions to forward (1) or reverse (-1) steps
cfg_def.tt = [0 -1 1 NaN; 1 0 NaN -1; -1 NaN 0 1; NaN 1 -1 0]; 
% example: given a transition from high-high (1) to high-low (2),
% cfg_def.tt(1, 2) is -1, so reverse step
  
cfg = ProcessConfig(cfg_def, cfg_in);

angle_tsd = state_tsd;
angle_tsd.data = NaN*angle_tsd.data;
angle_tsd.data(1) = 0; % assume we start at 0


% loop over state inputs; look up step given transition between current
% state and previous state at each step
for iState = 2:length(state_tsd.data)
   
    angle_tsd.data(iState) = cfg.tt(state_tsd.data(iState-1), state_tsd.data(iState));
    
end    % units of "angle" here are actual integer values that represent identical forward or backward increments of the wheel (i.e. not an actual angle value)


angle_tsd.data = cumsum(angle_tsd.data);
wheel_tsd.tvec = angle_tsd.tvec;  
wheel_tsd.data = angle_tsd.data./cfg.angle_to_turn; 