% Z = the sorted list of sessions by chronological order
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
Z = fd(Tnew);