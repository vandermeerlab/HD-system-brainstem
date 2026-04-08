function [m,stats] = typeI_II_singleSession_all_sessions(startSess, endSess, fd, cfg_in)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
cfg_def.startSess = startSess;
cfg_def.endSess = endSess;

cfg_out = ProcessConfig2(cfg_def, cfg_in);
if isempty(fd)
    fd = FindFiles('*keys.m');
end
if isempty(cfg_out.startSess); cfg_out.startSess = 1; end
if isempty(cfg_out.endSess); cfg_out.endSess = length(fd); end
for iSess = cfg_out.startSess : cfg_out.endSess
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN; disp(SSN);
    cfg_in = [];
    [m{iSess},stats{iSess},~] = typeI_II_singleSession(cfg_in);
    popdir;
end

