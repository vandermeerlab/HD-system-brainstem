function [AHVtsd] = GetAHV_values(orientationtouse, varargin)
%2020-01-14. JJS. This function uses dxdt to calculate velocity from orientationtouse data, akin to taking the derivative of position. 
%   Inputs:     orientationtouse - a TSD of orientationtouse values (should vary between -90 deg and 90 deg), calculated with GetOrientation.m
%   Outputs:    AHV         - Angular Head Velocity (for the headfixed mouse). Derivative of orientationtouse. Units of degrees per second.  

window = 0.1;
postsmoothing = .05;
extract_varargin; 

tic 
AHV = dxdt(orientationtouse.tvec, orientationtouse.data, 'window', window, 'postsmoothing', postsmoothing);
toc

AHVtsd = tsd(orientationtouse.tvec, AHV);

end

