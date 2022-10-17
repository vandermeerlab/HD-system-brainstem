function [STtstart, STtend, tlistX] = findSteadyPeriodsManual()
%2022-05-19. JJS. Manually identify start and end times for stationary periods in a recording session with platform veloicty info (AHV values).
%   code is borrowed from detectSaccadesManualCheck4.m

SSN = HD_GetSSN; disp(SSN);
%% check to see if eyetracking video exists and load pupil data
if exist(strcat(SSN, '-VT1.mp4')) ==2
    
    %% Load Events and get the session start time
    events_ts = LoadEvents([]);
    index = strfind(events_ts.label, 'Starting Recording');
    if index{1} == 1                                 % Start Recording should be in the first or second .label position.
        starttime = events_ts.t{1}(1);  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
    elseif index{2} == 1
        starttime = events_ts.t{2}(1); % for session with laser events, the start recording time moves to position 2.
    else
        error('could not find start time for this session')
    end
    %% Get timestamps from the .smi file
    [~, b, c] = fileparts(FindFile('*VT1.smi'));
    fn = strcat(b,c);
    tvec_raw = read_smi(fn);
    tvec = tvec_raw - starttime;
    
    %% Get the pupil trace dervied from Facemap
    f = FindFiles('*VT1_proc.mat');
    load(f{1}, 'pupil');
    pupilH = pupil{1}.com(:,2);
    %             pupilV = pupil{1}.com(:,1);
    
    % Subtract the mean
    meanH = nanmean(pupilH);
    %             meanV = nanmean(pupilV);
    % Make it into a TSD
    tsdH = tsd(tvec, pupilH - meanH);   % tsd of horizontal pupil position
    %             tsdV = tsd(tvec, pupilV - meanV);   % tsd of vertical pupil position
    
    diffH = tsd(tvec(2:end)', diff(tsdH.data)');     % Should this be (1:end-1) or (2:end)?
    %             diffV = tsd(tvec(2:end)', diff(tsdV.data)');     % tsd of vertical pupil velocity
    diffH.cfg.hdr{1}.Fs = 1 / median(diff(diffH.tvec));   % append the sampling rate
    
    %             tstart = diffH.tvec(1);
    %             tend = diffH.tvec(end);
    
    yyaxis left
    plot(diffH.tvec, diffH.data);
    ylabel('Horizontal Pupil Position')
    hold on
end

[AHV_tsd] = Get_AHV([]);

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
yyaxis right
plot(AHV_tsd.tvec, AHV_tsd.data); hold on
ylabel('AHV')

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

