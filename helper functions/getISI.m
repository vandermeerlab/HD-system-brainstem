function I = getISI(S, varargin)

% I = getISI(S)
%
% INPUT
% S is a ts of spiketimes
%
% OUTPUT
% I is a ctsd of interspike interval times

dt = 0.001;
events_ts = LoadEvents([]);
assert(strcmp(events_ts.label{1}, 'Starting Recording'))=1;
t0 = events_ts.t{1}(1);
ind = strcmp(events_ts.label, 'Stopping Recording');
t1 = events_ts.t(find(ind));

process_varargin(varargin);

if ~isa(S,'ts')
	error('S must be a ts');
end

T = t0:dt:t1;
spikeTimes = S;
spikeTimes0 = spikeTimes(1:(end-1));
spikeTimes1 = spikeTimes(2:end);
spikeDiffs = spikeTimes1 - spikeTimes0;
if any(spikeDiffs <= 2*dt)
	warning('ISIs exist < 2 dt.');
	I = ctsd(t0,dt,nan);
	return
end

spikeDiffs0 = [spikeDiffs; inf];
spikeDiffs1 = [inf; spikeDiffs];
spikeTimes2 = [spikeTimes+dt/2; spikeTimes-dt/2];
spikeDiffs2 = [spikeDiffs0; spikeDiffs1];

[spikeTimes2,order] = sort(spikeTimes2);
spikeDiffs2 = spikeDiffs2(order);

I = interp1(spikeTimes2', spikeDiffs2', T, 'nearest', 'extrap');
I = ctsd(t0,dt,I');
