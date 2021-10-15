function h = ShadedErrorbar(X, Y, E, varargin)

% h = ShadedErrorbar(X, Y, E, parms)
% 
% INPUTS 
%   X = x-data
%   Y = y-data
%   E = errors (assumed to be symmetric)
%   
% PARMS
%   color = 'k'
%   marker = '.'
%   L = E
%   U = E
%  
% If you want non-symmetric errorbars, redefine L and U.

color = 'k';
marker = '.';
L = E;
U = E;
process_varargin(varargin);

% make sure all are n x 1
if size(X,1)==1, X = X'; end
if size(Y,1)==1, Y = Y'; end
if size(L,1)==1, L = L'; end
if size(U,1)==1, U = U'; end

H = ishold; 
hold on

keep = ~isnan(X + L + U);
h0 = area(X(keep), [Y(keep)-L(keep) U(keep)+L(keep)]);
delete(h0(1));
set(h0(2), 'FaceColor', color, 'EdgeColor', 'None');
alpha(get(h0(2), 'children'), 0.25);
h(2) = h0(2);
h(1) = plot(X, Y, 'Color', color, 'Marker', marker, 'LineWidth', 2);

if H; hold on; else hold off; end
	


