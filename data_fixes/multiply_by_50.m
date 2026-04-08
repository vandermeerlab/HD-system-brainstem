function [fd_out, K, dt, samp] = multiply_by_50(fd, varargin)
% Resave the eye velocity variable in saccades-edited.mat so that it is correctly divided by dt (i.e. multiplied by a factor of 50).

if isempty(fd)
    fd = FindFiles('*keys.m');
end
disp(length(fd))
startSess = 1;
endSess = length(fd);
process_varargin(varargin);

for iSess = startSess: endSess
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN; disp(SSN);
    filename = strcat(SSN, '-saccades-edited.mat');
    if exist(strcat(SSN, '-saccades-edited.mat')) == 2
        load(strcat(SSN, '-saccades-edited.mat'));
        K(iSess) = k;
        disp(num2str(k))
        dt(iSess) = median(diff(diffH.tvec)); % this is the timestep (dt) for this session.
        samp(iSess) = 1/dt(iSess);
        diffHdeg = tsd(diffHdeg.tvec, diffHdeg.data*(samp(iSess)));  % 1/dt is the sampling rate for the camera for this session
        save(filename, 'diffHdeg', "-append")
        disp('data saved to saccades-edited.mat file')
    end
    fd_out{iSess} = fd{iSess};
end    
