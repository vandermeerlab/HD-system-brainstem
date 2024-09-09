function [csc_tsd, orientation, starttime, endtime, samplingrate, dt] = GetOrientationValues(cfg_in)
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
%           csc_tsd - a copy of the raw encoder CSC
%2020-03-04.  Updated to be compatible with Mvdm lab codeset.

cfg_def.rangetouse = 360;  % this was previously set to 180, which is incorrect. The full range of the platform is 180 degrees in either direction, which is 360 degrees total.
cfg_def.CheckPlot = 0;

dateswitch = datetime('2020-10-06');               % On this date I swtiched from using CSC21 to CSC33 for the platform encoder.
SSN = HD_GetSSN;
sessiondate = SSN(6:15);
sessiondate = datetime(sessiondate);

if sessiondate < dateswitch
    CSCtoUse = 21;     % CSC21 was used for the platform orientation encoder up until ______  % CSC33 was used after that date
else
    CSCtoUse = 33;
end

cfg = ProcessConfig2(cfg_def, cfg_in);

cfg_csc = [];
cfg_csc.fc = {FindFile(strcat('*CSC', num2str(CSCtoUse), '.ncs'))};
cfg_csc.VoltageConvFactor = 10^6;
csc_tsd = LoadCSC(cfg_csc);
starttime = csc_tsd.tvec(1);
endtime = csc_tsd.tvec(end);
csc_tsd.tvec = csc_tsd.tvec - csc_tsd.tvec(1); % MvdM: this probably shouldn't happen here because you need to know this information when loading spikes.

baseline = csc_tsd.data(1);               % It is critical that the rig is oriented in the same position from day to day when the arduino is turned on.
% Furthermore, the rig should be left to stand still for at least 30 sec or so to get a stable starting value to accurately subtract the baseline.
maxR = max(csc_tsd.data);                   % Clockwise rotation produces a less negative (positive) voltage change
maxL = min(csc_tsd.data);                   % Counterclockwise rotation produces a more negative (negative) voltage change
% ***Convention in Taube papers is for CW rotation to plotted as negative and CCW as positve, so I need a sign change
Rrange = abs(maxR-baseline);
Lrange = abs(maxL-baseline);
rangediff = Lrange/Rrange;
%                                          for rare sessions in which I do not move the platform through its full range, the 'fullrange' variable can be estimated as Fullrange =120550.
Fullrange = abs(maxL - maxR);           % ***Figure out what this value is and make sure it is the same from session to session. Add a warning.
format bank
disp(strcat('Fullrange = ', num2str(Fullrange)));
disp(strcat('Lrange = ', num2str(Lrange)));
disp(strcat('Rrange = ', num2str(Rrange)));
disp(strcat('Percentage difference equals__ ', num2str(100*(1-rangediff)), '%'));

% if Rrange ~= Lrange; warning('Left and Right ranges do not coincide.'); end

subtractedvoltage = tsd(csc_tsd.tvec, csc_tsd.data - baseline);
divisionconstant = Fullrange/cfg.rangetouse;  % this is the constant value to normalize by to get a max of 180 degrees rotation in either direction.
orientation = tsd(csc_tsd.tvec, subtractedvoltage.data./divisionconstant);
orientation.data = -orientation.data;  % THIS STEP IS CRITICAL. CW turns are defined as negative changes in angle (as with the unit circle). We need a sign change of the voltage values to make the raw voltage correspond to head direction (i.e. 'orientation').
dt = median(diff(orientation.tvec));
samplingrate = 1/dt;

if cfg.CheckPlot == 1
    figure
    subplot(2,1,1);
    plot(csc_tsd.tvec, csc_tsd.data);
    ylabel('Raw voltage')
    subplot(2,1,2);
    plot(orientation.tvec, orientation.data);
    ylabel('Orientation (deg)')
end
