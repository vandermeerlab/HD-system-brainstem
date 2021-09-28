function [PETH, outputS, outputT, outputI] = spike_PETHhd(S,t,varargin)

% spikePETH(S, t, varargin)
%
% input: 
%        TS S
%        Event times t
% 
% output:
%        TSD PETH - tsd with time as window, data as sum by binsize
% 
% if no outputs, plots
%
% parameters:
%   window = [-2 5]
%   dt = 0.001; % seconds
%   showHist = false
%   showSpikes = false
%
% ADR 2012-01
%
% ADR & JS 2012-11 fixed bug that histc changes its orientation when
% only one spike is in the window

%--------------------------
% parameters
%--------------------------
window = [-2 5];
dt = 0.01;
showSpikes = false;
showHist = false;
colortouse = 'b';
process_varargin(varargin);

%-------------------------
% prep
%--------------------------
nT = length(t);

outputS = [];
outputT = [];
outputI = [];

ISI = tsd(S.range(), [nan; diff(S.range())]);

%---------------------------
% go
%---------------------------

for iT = 1:nT
	S0 = data(  S.restrict(t(iT)+window(1), t(iT)+window(2)));
	I0 = data(ISI.restrict(t(iT)+window(1), t(iT)+window(2)));
	outputS = cat(1, outputS, S0-t(iT));
	outputT = cat(1, outputT, repmat(iT, size(S0)));
	outputI = cat(1, outputI, I0);
end

x = window(1):dt:window(2);
if length(outputS)==1
    H = histc(outputS, x)';
else
    H = histc(outputS, x);
end    
PETH = tsd(x', H);

if showSpikes
	% display
% 	figure	
	if ~isempty(outputS)
		plot(outputS, outputT,'.', 'color', colortouse);
	end
	axis([window(1) window(2) 0 nT]);
	xlabel('peri-event (sec)');
	ylabel('Event #');	
end

if showHist
	% display
	figure
	bar(PETH.range(), PETH.data());
	set(gca, 'XLim', [window(1) window(2)]);
	xlabel('peri-event (sec)');
	ylabel('Count');	
end