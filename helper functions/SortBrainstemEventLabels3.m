function [start_time, stop_time, laser_on, laser_off, arraysize, difflaser] = SortBrainstemEventLabels3

% 2023-07-21. JJS.
% Pull out Laser Start times and Laser Off times from the Events file.
% The actual labels for laser on/off TTL events varies. Therefore, it has to be handled with several different cases.

% Most sessions have laser events with the codes
%       'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0004)' &  'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0000)', or
%       'Laser On'                                                &  'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0000)'


events_ts = LoadEvents([]);
SSN = HD_GetSSN; disp(SSN); 

if length(events_ts.label) == 2    % These are sessions in which the only events are Acquisition START and STOP. No laser events.
    start_time = NaN;
    stop_time = NaN;
    laser_on = NaN;
    laser_off = NaN;
    arraysize = NaN;
    difflaser = NaN;
    disp('no laser events found')
    return
end

%% weird sessions with different TTL values
% these sessions are 'M094_2020_12_28', 'M096_2021_01_05', 'M112_2021_02_10'.
% Different values b/c of something I did with the master 8 settings. They have no opto cells. Will be ignored.
wrapper = @(events_ts) strcmp(events_ts, 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0008).');   % the few sessions with these values have no opto cells
A = cellfun(wrapper, events_ts.label);
bitEight_label = find(A); % index which label corresponds to (0x0008)

wrapper = @(events_ts) strcmp(events_ts, 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0010).');  % the few sessions with these values have no opto cells
A = cellfun(wrapper, events_ts.label);
bitTen_label = find(A); % index which label corresponds to (0x0008)

if ~isempty(bitEight_label) || ~isempty(bitTen_label)
    start_time = NaN;
    stop_time = NaN;
    laser_on = NaN;
    laser_off = NaN;
    arraysize = NaN;
    difflaser = NaN;
    disp('Session with different laser event labels: x0008 or x0010. Skipping session.')
    return
end

%% Grab the start time and stop times for the session. There may be one or two sessions with more than one start/stop time

wrapper = @(events_ts) strcmp(events_ts, 'Starting Recording');
A = cellfun(wrapper, events_ts.label);
start_label = find(A); % index which label corresponds to 'Start Recording'
start_time = events_ts.t{start_label};
% evt.label{1} = 'Starting Recording';
% evt.t{1} = start_time;

wrapper = @(events_ts) strcmp(events_ts, 'Stopping Recording');
A = cellfun(wrapper, events_ts.label);
stop_label = find(A); % index which label corresponds to 'Stopping Recording'
if ~isempty(stop_label)
    stop_time = events_ts.t{stop_label};
else
    warning('this session lacks a stop time. Could be a cheetah crash.')
end

%% Get event times
wrapper = @(events_ts) strcmp(events_ts, 'Laser On');
A = cellfun(wrapper, events_ts.label);
laser_on_label = find(A); % index which label corresponds to 'Laser On'
if ~isempty(laser_on_label)
    laser_on_unix = events_ts.t{laser_on_label};
    laser_on = laser_on_unix - start_time;
end

wrapper = @(events_ts) strcmp(events_ts, 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0000).');
A = cellfun(wrapper, events_ts.label);
bitZero_label = find(A); % index which label corresponds to (0x0000)
if ~isempty(bitZero_label)
    bitZero = events_ts.t{bitZero_label} - start_time;
end

wrapper = @(events_ts) strcmp(events_ts, 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0004).');
A = cellfun(wrapper, events_ts.label);
bitFour_label = find(A); % index which label corresponds to (0x0004)
if ~isempty(bitFour_label)
    bitFour = events_ts.t{bitFour_label} - start_time;
end
%% Exclude sessions without laser events, or the few sessions with laser event labels that are different (x0008, x0010) and I know don't have opto cells.

if isempty(laser_on_label) && isempty(bitFour_label)
    start_time = NaN;
    stop_time = NaN;
    laser_on = NaN;
    laser_off = NaN;
    arraysize = NaN;
    disp('no laser events detected for this session')
    return
end

%%
% CASE #1
% Sessions in which laser_on are present always have 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0000).'  as the laser off TTL label.
if exist('laser_on')
    laser_off = bitZero;
end
%%
% CASE #2
if exist('bitZero') && exist('bitFour')
    laser_on = bitFour;
    laser_off = bitZero;
end

%% Special cases
if strcmp(SSN, 'M281-2021-12-23') == 1    % very specific exception to this one session that crashed during a laser stim
    laser_on = laser_on(1:24);
    stop_time = NaN;
    arraysize = NaN;
end
if strcmp(SSN, 'M402-2022-11-02-1') == 1
    bitZero = bitZero(2:end);  % first value is whack
end

arraysize = min(length(laser_on), length(laser_off));
if length(laser_on) ~= length(laser_off)
    warning('unequal number of laser events')
end
for iEvent = 1:arraysize
    idiff(iEvent) = laser_off(iEvent) - laser_on(iEvent);
    try
        assert((idiff(iEvent)) < 1.01 &&  (idiff(iEvent) > 0.99)); % all laser intervals should be one second long. Check that they are indeed 1 second intervals.
    catch
        disp(iEvent)
        warning('laser duration is not correct')
    end
end

difflaser = laser_off(1:arraysize) - laser_on(1:arraysize); 


%% I'm not sure if I actually need this part
% if exist('laser_on')
%     laser_off = bitZero;
%     assert(sum(laser_on > bitZero) == 0)
%     bit0 = 1;
%     bit4 = 0;
% else
%     diff1 = bitZero(1:arraysize) - bitFour(1:arraysize);  % bitZero = laser_off; bitFour = laser_on
%     diff1mode = mode(diff1);
%     diff2 = bitFour(1:arraysize) - bitZero(1:arraysize);  % bitFour = laser_off; bitZero = laser_on
%     diff2mode = mode(diff2);
%
%     if diff1mode > 0 && diff2mode < 0
%         laser_off = bitZero;
%         laser_on = bitFour;
%         bit0 = bit0 + 1; % keep track of how many sessions are ordered this way, where 0x0000 = 'off'
%         assert(sum(laser_off > laser_on) == length(laser_on))   % the sum should be all ones [logical 'yes' values] and therefore equal to the length
%     elseif diff2 > 0 && diff1 < 0
%         laser_off = bitFour;
%         laser_on = bitZero;
%         bit4 = bit4 + 1; % keep track of how many sessions are ordered this way, where 0x0000 = 'on'
%         %         assert(sum(laser_off > laser_on) == length(laser_on))   % the sum should be all ones [logical 'yes' values] and therefore equal to the length
%     else
%         error('laser times do not make sense')
%     end
% end
%%