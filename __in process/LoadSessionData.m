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
%           sd.temporal - temporal saccade times
%           sd.nasal - nasal saccade times
%           sd.wheel -
%           sd.ExpKeys - structure with the elements loaded from LoadExpKeys. Should include start time [1x1], stop time [1x1], laser ON ts, Laser OFF ts, any manual entries, such as for DARK recording, optokinetic stim.
%           sd.cfg - contains config parameters that were used to generate the variables above
tic
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

SSN = HD_GetSSN; disp(SSN);
sd.SSN = SSN;
% -----------------------
% KEYS
% -----------------------
if Keys ==1
    keysfn = [strrep(SSN, '-', '_') '_keys'];
    assert(exist(keysfn, 'file')==2, 'Cannot find keys file %s.', keysfn);
    events_ts = LoadEvents([]);
    sd.ExpKeys = events_ts;
    sd.ExpKeys.SSN = SSN;
    sd.ExpKeys.fd = fd;
end
%-------------------------
% EVENTS
%-------------------------
if Events
    events_fn = fullfile(fd, [SSN '-Events.Nev']);
    assert(exist(events_fn, 'file')==2, 'Cannot find events file %s.', events_fn);
    events_ts = LoadEvents([]);
    sd.Events = events_ts;
end
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
    % New cheetah versions have timestamps that are in Unix Epoch Time
    index = strfind(sd.Events.label, 'Starting Recording');
    if index{1} == 1                                 % Start Recording should be in the first or second .label position.
        for iC = 1:length(S.t)
            S.t{iC} = S.t{iC} - events_ts.t{1}(1);  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
        end
    elseif index{2} == 1                   % this is for cases in which there is more than one start/stop time. Would not be usable for eye movement data, but can still be analyzed for AHV.
        for iC = 1:length(S.t)
            S.t{iC} = S.t{iC} - events_ts.t{2}(1);
        end
    else
        warning('could not find start record time')
    end
    sd.fn = {};
    for iC = 1:length(sd.fc)
        [~,sd.fn{iC}] = fileparts(sd.fc{iC});
        sd.fn{iC} = strrep(sd.fn{iC}, '_', '-');
    end
    sd.fn = sd.fn';
end
%-------------------------
% ANGULAR HEAD VELOCITY
%-------------------------
if AHV
    cfg = [];
    cfg.CheckPlot = 0;  % if CheckPlot == 1, will display the orientation plot.
    cfg.rangetouse = 360;  % this is the angular range for the platform, a critical value. Should be 360 for all sessions, except two, where it was set lower. *Note: hardcode thse here.
    [csc_tsd, orientation, samplingrate, dt] = GetOrientationValues(cfg); %#ok<ASGLU>
    subsample_factor = 10;
    orientationtousedata = downsample(orientation.data, subsample_factor);
    orientationtouserange = downsample(orientation.tvec, subsample_factor);
    
    window = 0.1; postsmoothing = .05;
    tic; AHV = dxdt(orientationtouserange, orientationtousedata, 'window', window, 'postsmoothing', postsmoothing); toc;
    AHV = -AHV; % THIS STEP IS NECESSARY BECAUSE dxdt GIVES VALUES THAT ARE CORRECT, BUT WITH A SIGN FLIP.
    sd.AHV = tsd(orientationtouserange, AHV);
    CheckPlot = 1;
    if CheckPlot ==1
        clf;
        plot(sd.AHV.tvec, sd.AHV.data);
        ylabel('AHV (deg/sec)')
        xlabel('Times (sec)')
    end
end
%-------------------------
% EYE MOVEMENTS
%-------------------------
if EYE
    cfg = [];
    cfg.threshAdj  = 4;
    cfg.threshH = 10;      % pixel threshold for CW rotation and Nasal Saccades.
    cfg.threshL = -10;     % pixel threshold for CCW rotation and Temporal Saccades.
    
    cfg.scalingfactor = 1;  % for shrinking the pupil trace so its the same height as diffH
    cfg.artifactThresh = 4;  % units of pixels sq.
    cfg.doPlotThresholds = 0;  % plot the eye velocity trace with thresholds and identified saccades
    cfg.doPlotEverything = 0;
    [sd.temporalSaccades, sd.nasalSaccades, ~, ~, ~, ~, ~, ~, ~] = processPupilData2(cfg);
end


if dirpushed
    popdir;
end
toc