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