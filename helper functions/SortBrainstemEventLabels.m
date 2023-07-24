function [start_time, stop_time, laser_on, laser_off, bit0, bit4, arraysize] = SortBrainstemEventLabels

% 2022-10-25. JJS.
% Get the correct timestamps for the different session events
% evt = struct;
events_ts = LoadEvents([]);
SSN = HD_GetSSN;

if length(events_ts.label) == 2
    start_time = NaN;  
    stop_time = NaN;   
    laser_on = NaN;   
    laser_off = NaN;   
    bit0 = NaN;   
    bit4 = NaN;   
    arraysize = NaN;  
    disp(SSN)
    disp('no laser events found')
    return
end

bit0 = 0;
bit4 = 0;

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

%% weird sessions with different TTL values
% these sessions are 'M094_2020_12_28', 'M096_2021_01_05', 'M112_2021_02_10'.
% Different values b/c of something I did with the master 8 settings. They have no opto cells. Will be ignored.
wrapper = @(events_ts) strcmp(events_ts, 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0008).');   % the few sessions with these values have no opto cells
A = cellfun(wrapper, events_ts.label);
bitEight_label = find(A); % index which label corresponds to (0x0008)
if ~isempty(bitEight_label)
    bitEight = events_ts.t{bitEight_label};
end
wrapper = @(events_ts) strcmp(events_ts, 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0010).');  % the few sessions with these values have no opto cells
A = cellfun(wrapper, events_ts.label);
bitTen_label = find(A); % index which label corresponds to (0x0008)
if ~isempty(bitTen_label)
    bitTen = events_ts.t{bitTen_label};
end

if ~exist('laser_on') && ~exist('bitZero')
    laser_on = NaN;
    laser_off = NaN;
    bit0 = NaN;
    bit4 = NaN;
    arraysize = NaN;
    return
end

%%
if exist('bitEight')==1 %#ok<*EXIST>
    laser_on = NaN;
    laser_off = NaN; %#ok<*NASGU>
    bit0 = NaN;
    bit4 = NaN;
    arraysize = NaN;
    return
end
if exist('bitTen')==1
    laser_on = NaN;
    laser_off = NaN;
    bit0 = NaN;
    bit4 = NaN;
    arraysize = NaN;
    return
end

if strcmp(SSN, 'M281-2021-12-23') == 1    % very specific exception to this one session that crashed during a laser stim
    laser_on = laser_on(1:24);
    stop_time = NaN;
    bit0 = NaN;
    bit4 = NaN;
    bitFour = 0;
    arraysize = NaN;
end
if strcmp(SSN, 'M402-2022-11-02-1') == 1
    bitZero = bitZero(2:end);  % first value is whack
end
    
% if strcmp(SSN, 'M112-2021-02-10') == 1
%     bitZero = bitZero(5:24);
%     bitFour = bitFour(5:24);
% end
    
   
Zerolength = length(bitZero);
    Fourlength = length(bitFour);
    arraysize = min([Zerolength Fourlength]);  % if the session ended before the last laser pulse turned off, then laser_off will be one short. 'arraysize' is designed to throw out the last laser pulse

%% I'm not sure if I actually need this part
if exist('laser_on')
    laser_off = bitZero;
    assert(sum(laser_on > bitZero) == 0)
    bit0 = 1;
    bit4 = 0;
else
    diff1 = bitZero(1:arraysize) - bitFour(1:arraysize);  % bitZero = laser_off; bitFour = laser_on
    diff1mode = mode(diff1);
    diff2 = bitFour(1:arraysize) - bitZero(1:arraysize);  % bitFour = laser_off; bitZero = laser_on
    diff2mode = mode(diff2);

    if diff1mode > 0 && diff2mode < 0
        laser_off = bitZero;
        laser_on = bitFour;
        bit0 = bit0 + 1; % keep track of how many sessions are ordered this way, where 0x0000 = 'off'
        assert(sum(laser_off > laser_on) == length(laser_on))   % the sum should be all ones [logical 'yes' values] and therefore equal to the length
    elseif diff2 > 0 && diff1 < 0
        laser_off = bitFour;
        laser_on = bitZero;
        bit4 = bit4 + 1; % keep track of how many sessions are ordered this way, where 0x0000 = 'on'
        %         assert(sum(laser_off > laser_on) == length(laser_on))   % the sum should be all ones [logical 'yes' values] and therefore equal to the length
    else
        error('laser times do not make sense')
    end
end
%

if exist('laser_on')
    laser_off = bitZero;
end

if length(laser_on) ~= length(laser_off)
    warning('unequal number of laser events')
end

laser_on = laser_on(1:arraysize);
laser_off = laser_off(1:arraysize);


    
    
    
    
    
    
    

