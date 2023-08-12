function [S] = LoadSpikesJeff(cfg)

% cfg = [];
cfg.uint = '64';
% spikefiles = FindFiles('*.t');
% cfg.fc = {spikefiles{iCell}};
S = LoadSpikes(cfg);

% New cheetah versions have timestamps that are in Unix Epoch Time
events_ts = LoadEvents([]);

wrapper = @(events_ts) strcmp(events_ts, 'Starting Recording');
A = cellfun(wrapper, events_ts.label);
Startindex = find(A); % index which label says 'Start Recording'
starttime = events_ts.t{Startindex}(1); % use the very first start record time

for iC = 1:length(S.t)
    S.t{iC} = S.t{iC} - starttime;  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
end

