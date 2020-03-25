function [S, Q, F] = plotRawTuningCurve(varargin)
%2020-03-17. JJS. Load spikes, calculate firing rate, and plot a 'raw version' of Firing Rate by AHV values. 
%   Detailed explanation goes here

extract_varargin;

[csc_tsd, orientation, ~, ~] = GetOrientationValues([]); %#ok<ASGLU>
[orientationtouse] = downsampleOrientationValues(orientation, 10);
[AHVtsd] = GetAHV_values(orientationtouse);
% newrange = AHVtsd.tvec - AHVtsd.tvec(1);    % timestamps are in microseconds instead of seconds. Figure out why this is still happenning.
% AHVtsdnew = tsd(newrange, AHVtsd.data);
S = LoadSpikes([]);
% tc_out = TuningCurves([], S, AHVtsd);



end

