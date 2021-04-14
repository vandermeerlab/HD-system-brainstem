function [S] = LoadSpikesJeff

cfg = [];
cfg.uint = '64';
% spikefiles = FindFiles('*.t');
% cfg.fc = {spikefiles{iCell}};
S = LoadSpikes(cfg);

% New cheetah versions have timestamps that are in Unix Epoch Time
events_ts = LoadEvents([]);
assert(strcmp(events_ts.label{1}, 'Starting Recording'))=1;
for iC = 1:length(S.t)
    S.t{iC} = S.t{iC} - events_ts.t{1}(1);  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
end
