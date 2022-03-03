function [d, speed, cfg] = ConvertWheeltoSpeed(cfg_in, wheel_tsd)

% JJS. 2020-03-03. 
% Take the wheel revolution data (number of wheel revolutions as a function of time) and create a tsd of distance travelled and instantaneous speed. 

cfg_def = [];
cfg_def.downsample = 20; 
cfg_def.wheel_diam = 15; % diameter of the wheel in centimeters. Measured on 3/3/2022. 
  
cfg = ProcessConfig(cfg_def, cfg_in);

distance = wheel_tsd.data * pi * cfg.wheel_diam; 

timestamps = downsample(wheel_tsd.tvec, cfg.downsample); 
distance = downsample(distance, cfg.downsample); 

% distance = cat(2, 0, distance);  % add a zero value so that when we diff, the vector length remains the same as before. 

d.tvec = timestamps;
d.data = distance; 

speed.tvec = timestamps;
speed.data = dxdt(d.tvec, d.data);

end