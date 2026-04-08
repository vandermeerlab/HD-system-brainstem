function [STtstart, STtend, tlistX] = findSteadyPeriods()
%2022-05-19. JJS. Manually identify start and end times for stationary periods in a recording session with platform veloicty info (AHV values).
%   code is borrowed from detectSaccadesManualCheck4.m 

[AHV_tsd] = Get_AHV([]);

SSN = HD_GetSSN; disp(SSN);
if exist(strcat(SSN, '-AHV_StationaryTimes.mat'), 'file') == 2
    disp('-AHV_StationaryTimes.mat already exists');
else
    disp('no StationaryTimes file present')
end
f = input('Do you want to overwrite (o) the existing file/create a new file, append (a), or skip (s) this session? [o/a/s]', 's');
if strcmp(f, 'o')
    disp('File will be overwritten')
elseif strcmp(f, 'a')
    overwrite = 0;
    disp('Selections will be appended to the existing data')
elseif strcmp(f, 's')
    skip = 1;
else
    derror('unrecognized response')
end
clf
plot(AHV_tsd.tvec, AHV_tsd.data); hold on

%% While LOOP
fprintf(1, '\n');
fprintf(1, '\n');
disp('start and stop times for STATIONARY periods')
count = 0;
tlistX= [];
tlistY = [];

while(1)
    fprintf(1, '\n');
    feedback = input('Do you want to continue, y/n? ...','s');
    if feedback == 'n'
        break
    end
    clear feedback
    count = count + 1;
    fprintf(1, '\n');
    disp('Zoom into region of interest. Press return to stop zoom.')
    zoom on;
    pause() % you can zoom with your mouse and when your image is okay, you press any key
    zoom off; % to escape the zoom mode
    fprintf(1, '\n');
    disp('Select a PAIR of points to indicate a stationary EPOCH.')
    [x,y, ~] =ginput(2);  % Choose two points to mark the START/END of a stationary period. 
    plot(x(1), y(1), 'g.', 'MarkerSize', 25)
    plot(x(2), y(2), 'r.', 'MarkerSize', 25) 
    accept = input('Accept these selections?  enter/n ...', 's');
    if strcmp(accept,'') || isempty(accept)
        tlistX(end+1:end+length(x)) = x;
        tlistY(end+1:end+length(y)) = y;
        disp('points selected')
    else
        disp('Continue through prompts to re-select points')
        plot(x, y, 'w.', 'MarkerSize', 25)
    end
end

assert(mod(length(tlistX),2) == 0); 

STtstart = tlistX(1:2:end);
STtend = tlistX(2:2:end);

disp('here are the start and end times')
STtstart %#ok<*NOPRT>
STtend

save(strcat(SSN, '-AHV_StationaryTimes.mat'), 'STtstart', 'STtend', 'tlistX');
disp('data saved') 


end

