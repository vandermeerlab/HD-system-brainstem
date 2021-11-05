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

Keys = true; %1 = load keys.m, 0 = don't
VT1 = true;  %1 = load VT1, 0 = don't
VT2 = true;  %1 = load VT2, 0 = don't
Spikes = true;  %1 = load spikes, 0 = don't
% HF  = true;  %1 = load *DD.mat, 0 = don't.  .mat file doesn't exist yet.
Use__Ts = false; % load ._t cells
Events = true; % load events file
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


















if dirpushed
    popdir;
end