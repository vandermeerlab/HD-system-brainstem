function A=scatterplotc(x,y,C,varargin)
%  H=scatterplotC(x,y,C)
%  H=scatterplotC(x,[],C)
%
% Plots x and/or y with color proportional to C. 
%
%INPUTS
%  x,y -- input vector(s) to plot
%  C   -- color level to plot         
%
%OUTPUTS
%  H   -- figure handle
%
%Parameters\Defaults
%
%  plotchar   = 'o'
%  NumColors  = 100;
%  solid_face = false;
%  crange = [min(C) max(C)]; % caxis
%  displayAll = if true C(C>maxC) shows up as maxC, C(C<minC) as minC, 
%               else points outside the range are not displayed
% 
% if you want log-scale, use scatterplotc(x,y,log(C));
% 
%
% JCJ 25 Oct 2002
% JCJ 01 May 2005 -- enabled varargin for use of all line graphic properties as in plot
%                 -- also included ability to make solid markers
% ADR updated for current system
% ADR 3 Nov 2012 fixed outofrange error

%-------------------
% DEFAULTS
%-------------------
plotchar   = 'o';
NumColors  = 100;
solid_face       = false;
crange = [];
displayAll = true;

varargin = process_varargin(varargin);

%-------------------
% CHECK
%-------------------
assert(min(size(x))==1, 'x must be vector.');
assert(isempty(y) || min(size(y))==1, 'y must be vector.');

if isempty(y) % called as (x,[],C)
	y = x;
	x = 1:length(y);
end

assert(all(size(x)==size(y)), 'x,y,C must be same size.');
assert(all(size(C)==size(x)), 'x,y,C must be same size.');

%--------------------
% build color matrix
%--------------------
if isempty(crange)
	crange = [min(C) max(C)];
end

% displayall
if displayAll
C(C>crange(2)) = crange(2);
C(C<crange(1)) = crange(1);
end

indexC = 1 + floor(((C-crange(1))/(crange(2)-crange(1)))*(NumColors-1));
    
caxis([0 NumColors]);
cmap = colormap(jet(NumColors));          % Make 'Hot' colormap for FR coding

%--------------------
% get figure
%--------------------
F    = gcf;                               % Gets figure handle to plot on; if no figure, makes one
if ishold()
	holdMemory = true; hold on;
else
	clf; holdMemory = false; hold on;
end

%-------------------
% plot
%--------------------
for iColor = 1:NumColors
	FR2plot = (indexC==iColor); % Find array positions of points
	if sum(FR2plot)
		h = plot(x(FR2plot),y(FR2plot),plotchar,'Color', cmap(iColor,:),varargin{:});
		if solid_face         %                 -- also included ability to make solid markers
			set(h, 'MarkerFaceColor', cmap(iColor,:));
		end
	end
end

%------------------
% unpack
%------------------
if holdMemory, hold on; else hold off; end

if nargout
    A=F;
end

