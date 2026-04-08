function [outputArg1,outputArg2] = calc_ipsi_contra(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;

for iSess = 1:numSess
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN; disp(SSN)
    if ~exist(strcat(SSN,'-VT1.mp4')); error('VIDEO FILE NOT FOUND'); end
    EvalKeys;
end

