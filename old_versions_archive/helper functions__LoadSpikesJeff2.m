function [S] = LoadSpikesJeff2
SSN = HD_GetSSN; 

cfg = [];
cfg.uint = '64';
% spikefiles = FindFiles('*.t');
% cfg.fc = {spikefiles{iCell}};
Sold = LoadSpikes(cfg);

% New cheetah versions have timestamps that are in Unix Epoch Time
events_ts = LoadEvents([]);

events_ts = LoadEvents([]);
wrapper = @(events_ts) strcmp(events_ts, 'Starting Recording');
A = cellfun(wrapper, events_ts.label);
Startindex = find(A); % index which label says 'Start Recording'
starttime = events_ts.t{Startindex}(1); % use the very first start record time
wrapper2 = @(events_ts) strcmp(events_ts, 'Stopping Recording');
B = cellfun(wrapper2, events_ts.label);
EndIndex = find(B);
endtime = events_ts.t{EndIndex};
% assert(strcmp(events_ts.label{1}, 'Starting Recording')==1)
for iC = 1:length(Sold.t)
    if strcmp(SSN, 'M054-2020-03-03') == 1
        S.t{iC} = Sold.t{iC} - Sold.t{iC}(1);
    else
        S.t{iC} = Sold.t{iC} - starttime;  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
    end
end


