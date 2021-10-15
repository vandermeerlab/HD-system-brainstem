function [S] = LoadSpikesJeff

cfg = [];
cfg.uint = '64';
% spikefiles = FindFiles('*.t');
% cfg.fc = {spikefiles{iCell}};
S = LoadSpikes(cfg);

% New cheetah versions have timestamps that are in Unix Epoch Time
events_ts = LoadEvents([]);
index = strfind(events_ts.label, 'Starting Recording');
if index{1} == 1;
    for iC = 1:length(S.t)
        S.t{iC} = S.t{iC} - events_ts.t{1}(1);  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
    end
elseif index{2} == 1;
    for iC = 1:length(S.t)
        S.t{iC} = S.t{iC} - events_ts.t{2}(1);  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
    end
else
    warning('could not find start record time')
end

