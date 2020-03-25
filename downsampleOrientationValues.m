function [orientationtouse] = downsampleOrientationValues(orientation, downsamplefactor)
%JJS. 2020-02-13.  Downsamples the tsd 'orientation' into a more maneagable size for use in the function dxdt. 
% 2020-03-04. Make compatible with vdM lab codebase. 

orientationtousedata = downsample(orientation.data, downsamplefactor);
orientationtouserange = downsample(orientation.tvec, downsamplefactor);

orientationtouse = tsd(orientationtouserange, orientationtousedata);


% change this to decimate 
end
