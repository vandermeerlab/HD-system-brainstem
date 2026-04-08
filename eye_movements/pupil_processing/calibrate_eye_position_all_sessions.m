function [k] = calibrate_eye_position_all_sessions(fd, cfg_in)
% JJS. 2024-07-25.
%   This function runs calibrate_eye_position_single_session_ver2 on all sessions, one after another, in the parent directory. See that function
% for a description of the procedure.

cfg_def.startSess = 1;
cfg_def.endSess = length(fd);
cfg_out = ProcessConfig(cfg_def,cfg_in);

if isempty(fd)
    fd = FindFiles('*keys.m');
end
for iSess = cfg_out.startSess : cfg_out.endSess
    pushdir(fileparts(fd{iSess}));
    %% Make sure that eye tracking exists for this session
    SSN = HD_GetSSN;
    if exist(strcat(SSN, '_pupil velocity calibrated.fig')) == 2  % if the fig exists, then this session has already been done. *Hacky. Change to see if variable k exists.
        disp('pupil calibration variables already exist. Skipping session')
    elseif exist(strcat(SSN, '-VT1.smi')) ~= 2
        disp('eye tracking data does not exist for this session. Skipping session.')
        k(iSess) = NaN;
    else
        sd = LoadSessionData([]);
        [k(iSess), ~, ~] = calibrate_eye_position_single_session_ver2(sd, cfg_out);
    end
    popdir;
end
