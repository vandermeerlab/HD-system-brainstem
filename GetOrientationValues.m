function [csc_tsd, orientation, samplingrate, dt] = GetOrientationValues(cfg, varargin)
%2019-12-20. JJS.
%   Pulls out the raw trace from the encoder and converts the signal into orientation in degrees.
%   Current configuration of the arduino is for it to process 180 deg rotation total (90 deg clockwise, 90 deg counterclockwise).

% Update: the arduino was reconfigured to process 360 deg of rotation total
%  
%
%   When the arduino is turned on, the angle of the rig at that time will correspond to a voltage directly in the middle of the range
%   of the signal (62.5mV). -180deg is ~0mV and +180 deg is ~125mV. The signal in cheetah roughly ranges between 0 and -20,000. This will depend
%   upon the gain of that channel in cheetah. It is likely? recorded in units of microvolts.

%  Inputs:  name (in quotes) of the CSC file that contains the encoder signal. Currently, the encoder is plugged in to AD Channel 0,
%           which corresponds to CSC21 in cheetah (bc of the mapping between the probe and omnetics connector).

%  Outputs: a tsd of orientation valuees with size equal to that of the original CSC. The range of values depends upon the setting of the arduino, and it needs
%  to be entered here as the variable 'rangetouse'.

%2020-03-04.  Updated to be compatible with Mvdm lab codeset.

rangetouse = 360;  % this was previously set to 180, which is incorrect. The full range of the platform is 180 degrees in either direction, which is 360 degrees total. 
CheckPlot = 0;
extract_varargin;

if isempty(cfg) == 1;
    cfg = [];
    cfg.fc = {FindFile('*CSC21.ncs')};
    cfg.VoltageConvFactor = 10^6;
end

csc_tsd = LoadCSC(cfg);
csc_tsd.tvec = csc_tsd.tvec - csc_tsd.tvec(1);

baseline = csc_tsd.data(1);               % It is critical that the rig is oriented in the same position from day to day when the arduino is turned on.
% Furthermore, the rig should be left to stand still for at least 30 sec or so to get a stable starting value to accurately subtract the baseline.
maxR = max(csc_tsd.data);                   % Clockwise rotation produces a less negative (positive) voltage change
maxL = min(csc_tsd.data);                   % Counterclockwise rotation produces a more negative (negative) voltage change
% ***Convention in Taube papers is for CW rotation to plotted as negative and CCW as positve, so I need a sign change
Rrange = abs(maxR-baseline);
Lrange = abs(maxL-baseline);

Fullrange = abs(maxL - maxR);           % ***Figure out what this value is and make sure it is the same from session to session. Add a warning.
disp(strcat('Fullrange = ', num2str(Fullrange)))
disp(strcat('Lrange = ', num2str(Lrange)))
disp(strcat('Rrange = ', num2str(Rrange)))

if Rrange ~= Lrange; warning('Left and Right ranges do not coincide.'); end

subtractedvoltage = tsd(csc_tsd.tvec, csc_tsd.data - baseline);
divisionconstant = Fullrange/rangetouse;  % this is the constant value to normalize by to get a max of 180 degrees rotation in either direction.
orientation = tsd(csc_tsd.tvec, subtractedvoltage.data./divisionconstant);
dt = median(diff(orientation.tvec));
samplingrate = 1/dt;

if CheckPlot == 1; 
    figure
    subplot(2,1,1);
    plot(csc_tsd.tvec, csc_tsd.data);
    subplot(2,1,2);
    plot(orientation.tvec, orientation.data);
    ylabel('Orientation (deg)')
end
