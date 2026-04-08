function [endtimetouse, csc_tsd] = findSessionEndTime
%2024-11-22. JJS.

dateswitch = datetime('2020-10-06');               % On this date I swtiched from using CSC21 to CSC33 for the platform encoder.
SSN = HD_GetSSN;
sessiondate = SSN(6:15);
sessiondate = datetime(sessiondate);
if sessiondate < dateswitch
    CSCtoUse = 21;     % CSC21 was used for the platform orientation encoder up until ______  % CSC33 was used after that date
else
    CSCtoUse = 33;
end

cfg_csc = [];
cfg_csc.fc = {FindFile(strcat('*CSC', num2str(CSCtoUse), '.ncs'))};
cfg_csc.VoltageConvFactor = 10^6;
csc_tsd = LoadCSC(cfg_csc);

starttime = csc_tsd.tvec(1); % subtraction is necessary b/c cheetah stores timestamps in Linux time.
endtime = csc_tsd.tvec(end);
endtimetouse = endtime - starttime;

csc_tsd = tsd(csc_tsd.tvec - starttime, csc_tsd.data);

