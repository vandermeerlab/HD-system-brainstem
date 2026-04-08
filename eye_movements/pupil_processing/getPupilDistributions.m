function [K] = getPupilDistributions(fd,cfg_in)
%2024-07-26. JJS. 
%   This function collects all of the k values (conversion factors for pixel to degrees) as well as the range of pupil positions and velocities for each session in the parent directory. 
cfg_def.startSess = 1;
cfg_def.endSess = length(fd);
cfg_out = ProcessConfig(cfg_def,cfg_in);

if isempty(fd)
    fd = FindFiles('*keys.m');
end
for iSess = cfg_out.startSess : cfg_out.endSess
    pushdir(fileparts(fd{iSess}));
%     SSN = HD_GetSSN; disp(SSN);
    clear k
    clear tsdHdeg
    clear diffHdeg
    load(strcat(SSN, '-saccades-edited.mat'));
    if exist('k', 'var') ==1 
        K(iSess) = k;
    else
        K(iSess) = NaN;
    end
    popdir
end

