function sd = LoadSessionData(fd, varargin)
% 2021-11-02. JJS.
% Loads all or most of the relevant information in a single session of head-fixed mouse brainstem recording. Must be in session folder to begin.
% Timestamps for all elements are centered at zero. This involves subtraction, because cheetah now uses unix time.
% For now, saccade times are calculated on the fly, although this is time intensive. May make sense to save these values once saccade code is close to optimal.

%  Inputs:
%           varargin - config variable.
%

%  Outputs:
%           sd.S - the init structure with all initialized variables
%           sd.AHV - tsd with AHV values (deg./sec) by time.
%           sd.temporal - temporal saccadegit a times
%           sd.nasal - nasal saccade times
%           sd.wheel -
%           sd.ExpKeys - structure with the elements loaded from LoadExpKeys. Should include start time [1x1], stop time [1x1], laser ON ts, Laser OFF ts, any manual entries, such as for DARK recording, optokinetic stim.
%           sd.cfg - contains config parameters that were used to generate the variables above

% subtractStartTime = 1; % New cheetah versions have timestamps that are in Unix Epoch Time. Not sure if some of my oldest sessions are on an earlier cheetah version where this is not necessary (i.e. starts at time = 0)
CheckOrientation = 0;
CheckAHV = 0;
Keys = true; %1 = load keys.m, 0 = don't
Events = true; % load events file
Spikes = true;  %1 = load spikes, 0 = don't
Use__Ts = false; % load ._t cells
AHV = true;
EYE = true;
process_varargin(varargin);

if ~isempty(fd)
    pushdir(fd);
    dirpushed = true;
else
    fd = pwd;
    dirpushed = false;
end
assert(exist(fd, 'dir')==7, 'Cannot find directory %s.', fd);

SSN = HD_GetSSN;
sd.SSN = SSN;
% -----------------------
% KEYS
% -----------------------
if Keys ==1
    keysfn = [strrep(SSN, '-', '_') '_keys'];
    if exist(keysfn, 'file')~= 2
        warning('Cannot find keys file %s.', keysfn);
    end
    EvalKeys(fd);
    sd.ExpKeys = ExpKeys;
end
%-------------------------
% EVENTS
%-------------------------
if Events
    events_fn = fullfile(fd, [SSN '-Events.Nev']);
    if exist(events_fn, 'file') ~= 2
        warning('Cannot find events file %s.', events_fn);
    end
    events_ts = LoadEvents([]);
    sd.Events = events_ts;
end
% New cheetah versions have timestamps that are in Unix Epoch Time, so you have to subtract the start time.

% wrapperA = @(events_ts) strcmp(events_ts, 'Starting Recording');
% A = cellfun(wrapperA, events_ts.label);
% wrapperB = @(events_ts) strcmp(events_ts, 'Stopping Recording');
% B = cellfun(wrapperB, events_ts.label);
% Startindex = find(A); % index which label says 'Start Recording'
% starttime = events_ts.t{1, Startindex}; % use the very first start record time
% Endindex = find(B);
% if strcmp(SSN, 'M281-2021-12-23') ==1   % This session crashed and doesn't have a cheetah endtime. JJS. 2024-07-26
%     endtime = 396.75;
% else
% endtime = events_ts.t{1, Endindex};
% end
% sd.starttime = starttime;
% sd.endtime = endtime;
% sd.SessLength = sd.endtime - sd.starttime;
%-------------------------
% ANGULAR HEAD VELOCITY
%-------------------------
if AHV
    cfg = [];
    cfg.CheckPlot = 0;  % if CheckPlot == 1, will display the orientation plot.
    cfg.rangetouse = 360;  % this is the angular range for the platform, a critical value.
    [csc_tsd, orientation, samplingrate, dt] = GetOrientationValues(cfg); %#ok<ASGLU>
    subsample_factor = 10;
    orientationtousedata = downsample(orientation.data, subsample_factor);
    orientationtouserange = downsample(orientation.tvec, subsample_factor);
    orientationtouse = tsd(orientationtouserange, orientationtousedata);
    sd.orientation = orientationtouse;
    sd.orientationsamplingrate = samplingrate/subsample_factor;
    
    if CheckOrientation ==1
        figure;
        plot(orientationtouserange, orientationtousedata);
        xlabel('Time (sec)')
        ylabel('Heading Direction (deg)')
    end
    window = 0.1; postsmoothing = .05;
    AHV = dxdt(orientationtouserange, orientationtousedata, 'window', window, 'postsmoothing', postsmoothing);
    AHV = -AHV; % THIS STEP IS NECESSARY BECAUSE dxdt GIVES VALUES THAT ARE CORRECT, BUT WITH A SIGN FLIP.
    sd.AHV = tsd(orientationtouserange, AHV);
    sd.AHV.dt = median(diff(sd.AHV.tvec));
    if CheckAHV ==1
        figure;
        plot(sd.AHV.tvec, sd.AHV.data);
        ylabel('AHV (deg/sec)')
        xlabel('Time (sec)')
    end
end
% Start and end session times
sd.starttime = csc_tsd.tvec(1); % uses the first and last timestamps of the encoder signal instead of the Cheetah start/stop times in the Events file (in case the session crashed)
sd.endtime = csc_tsd.tvec(end);
sd.SessLength = sd.endtime - sd.starttime;
%-------------------------
% SPIKES
%-------------------------
if Spikes ==1
    cfg = [];
    cfg.uint = '64';
    if Use__Ts == 1
        cfg.fc = cat(1, fc, FindFiles('*._t', 'CheckSubdirs',0));
    else
        cfg.fc = FindFiles('*.t', 'CheckSubdirs', 0);
    end
    sd.fc = cfg.fc;
    S = LoadSpikes(cfg);
    % Start Recording should be in the first or second .label position.
    for iC = 1:length(S.t)
        S.t{iC} = S.t{iC} - sd.starttime;  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
    end
    sd.S = S;
    
    sd.fn = {};
    for iC = 1:length(sd.fc)
        [~,sd.fn{iC}] = fileparts(sd.fc{iC});
        sd.fn{iC} = strrep(sd.fn{iC}, '_', '-');
    end
    sd.fn = sd.fn';
end
%-------------------------
% EYE MOVEMENTS
%-------------------------
if EYE
    success = 0; %#ok<NASGU>
    try
        load(FindFile('*saccades-edited.mat'), 'temporalSaccades', 'temporalAmplitudes', 'nasalSaccades', 'nasalAmplitudes', 'tsdH', 'tsdV', 'diffH', 'diffV')
        % ____saccades = timestamps.
        % ____amplitudes = saccade amplitdues.
        % tsdH is horizontal pupil position
        % tsdV is vertical pupil position
        success = 1;
        if success
            sd.temporalSaccades = temporalSaccades;
            sd.temporalAmplitudes = temporalAmplitudes;
            sd.nasalSaccades = nasalSaccades;
            sd.nasalAmplitudes = nasalAmplitudes;
            sd.tsdH = tsdH;
            sd.tsdH.data = sd.tsdH.data'; % TuningCurves() requires this format whereby tvec and data have opposite orientations (1 x n and n x 1, instead of 1 x n and 1 x n). 2023-07-19.
            sd.tsdH_dt = median(diff(tsdH.tvec));
            sd.tsdV = tsdV;
            sd.diffH = diffH;
            sd.diffV = diffV;
        end
    catch
        warning('cannot find saccade data')
    end
    % Separate MOVING vs. SPONTANEOUS saccades
    [~, ~, ~, ~, sd.nasal_timestamps_MOVING, sd.temporal_timestamps_MOVING] = isolateManualSaccades();
    [~, ~, sd.nasal_timestamps_REST, sd.temporal_timestamps_REST] = isolateStationarySaccades();
    
    sd.numNasal_moving = length(sd.nasal_timestamps_MOVING); if sd.numNasal_moving == 1; sd.numNasal_moving = 0; end  % if this array is a single NaN, then there are in fact no saccades.
    sd.numNasal_stationary = length(sd.nasal_timestamps_REST); if sd.numNasal_stationary == 1; sd.numNasal_stationary = 0; end % if this array is a single NaN, then there are in fact no saccades.
    sd.numTemporal_moving = length(sd.temporal_timestamps_MOVING); if sd.numTemporal_moving == 1; sd.numTemporal_moving = 0; end % if this array is a single NaN, then there are in fact no saccades.
    sd.numTemporal_stationary = length(sd.temporal_timestamps_REST); if sd.numTemporal_stationary == 1; sd.numTemporal_stationary = 0; end % if this array is a single NaN, then there are in fact no saccades.
end

%----------------------------
% WHEEL ENCODER (Quadrature)
%----------------------------
encoderCSC = strcat(SSN, '-CSC34.Ncs');
if exist(encoderCSC)
    updownTSD = getQEupdown([]);
    state_tsd = ConvertQEUpDownToState(updownTSD);
    [angle_tsd, wheel_tsd] = ConvertQEStatesToAngle([], state_tsd); %#ok<*ASGLU>
    [d, speed, cfg] = ConvertWheeltoSpeed([], wheel_tsd);  % speed the wheel (in centimeters per second)
    
    d.tvec = d.tvec - sd.starttime;
    speed.data = -speed.data;           % we want forward motion to be displayed as a positive velocity. This requires a sign reversal.
    speed.tvec = speed.tvec - sd.starttime;
    speed.dt = median(diff(speed.tvec));
    
    sd.d = d; % total distance travelled on the wheel (in centimenters)
    sd.speed = speed;
end
%----------------------------
% PLATFORM ROTATION PERIODS
%----------------------------
try
    success = 0; %#ok<*NASGU>
    filename = strcat(SSN, '-AHV_StationaryTimes.mat');
    load(filename);
    success = 1;
    sd.STstart = STtstart;
    sd.STend = STtend;
catch
    warning('cannot find platform stationary times file')
end

if dirpushed
    popdir;
end


%% alternate way of find start and stop recording times
% start_index = strfind(sd.Events.label, 'Starting Recording');
% if start_index{1} == 1                                 % Start Recording should be in the first or second .label position.
%     for iC = 1:length(S.t)
%         S.t{iC} = S.t{iC} - events_ts.t{1}(1);  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
%     end
% elseif start_index{2} == 1                   % this is for cases in which there is more than one start/stop time. Would not be usable for eye movement data, but can still be analyzed for AHV.
%     for iC = 1:length(S.t)
%         S.t{iC} = S.t{iC} - events_ts.t{2}(1);
%     end
% else
%     warning('could not find start record time')
% end
%%
% A = cell2mat(cellfun(@isempty, start_index,'UniformOutput',0));
% B = A == 0;
% indextouse = find(B);
% start_time = sd.Events.t{indextouse};
