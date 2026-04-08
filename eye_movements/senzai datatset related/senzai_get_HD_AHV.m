function [heading_degreesTSD, AHVtsd] = senzai_get_HD_AHV(mouseID, varargin)
%2024-10-31. JJS. Load the heading data, unwraps it, and calculates AHV. Also generates AHV TCs. 
%   Inputs: 
%               mouseID             - string with mouse identity
%               HD_neuronsToUse     - double. List of neurons with strong HD tuning. If empty, will skip calculating AHV TCs. 
%   Outputs:
%               heading_degreesTSD  - TSD of tvec and heading (in degrees of angle)
%               AHVtsd              - TSD of tvec and AHV (in degrees of angle per second)

% doPlot = 1;
process_varargin(varargin)
%% Get heading
filename = strcat('YutaTest', mouseID, '_OpenField_HeadDirection.mat');
heading = importdata(filename); % contains Neck&NoseOmitIdx [value of 1 = omit], rho(distance between nose to neck),theta (heading, in radians), t (timestamp, IN SECONDS)
heading_degrees = heading.theta*180/pi; % convert from radians to degrees
heading_degreesTSD = tsd(heading.t, heading_degrees');

%% Get AHV
tracking_dur = (heading.t(end) - heading.t(1))/60; % How many minutes long the tracking was in the open field. This is 98 minutes for M77.
disp(strcat('Session length was =', num2str(round(tracking_dur)), 'minutes'))
Q = unwrap(heading.theta);
Q_deg = Q*180/pi;
window = 0.1;
postsmoothing = .05;
AHV = dxdt(heading.t, Q_deg, 'window', window, 'postsmoothing', postsmoothing);
AHVtsd = tsd(heading.t, AHV);
% heading_sampling_rate = 1/median(diff(heading.t));  % should be 50Hz

end

