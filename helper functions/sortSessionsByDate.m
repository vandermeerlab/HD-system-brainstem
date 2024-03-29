function [fdSort, SSNonly] = sortSessionsByDate(fd) 
% Input:            fd 
% Output:           fdSort = the sorted list of sessions by chronological order
% 
datetouse = [];
for iF = 1:length(fd)
    temp = fd{iF}(33:42);
    datetouse = strvcat(datetouse, temp); 
    
end
% datetimetouse{iF} = datetime(datetouse{iF});
time = datetime(datetouse);
x = 1:length(datetouse); x = x';
t = table(time, x);
T = sortrows(t, 'time'); 
Tnew = table2array(T(:,2));
fdSort = fd(Tnew);

for iF = 1:length(fdSort)
    [~, b, ~] = fileparts(fdSort{iF});
    SSNonly{iF} = b; 
end
SSNonly = SSNonly';
