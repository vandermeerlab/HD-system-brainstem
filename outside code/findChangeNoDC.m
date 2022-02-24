function inds = findChangeNoDC(varargin)
%findChangeNoDC function  inds = findChangeNoDC(vec,threshold,negation(optional),plot(optional))- returns the points where there was 
% a change in input vector which is bigger than threshold.
%   inputs: 
%           1. the 1D vector.
%           2. the threshold
%           3. include negative responses if = 2;
%           4. binary value if one add plots

vec = varargin{1};
thresh = varargin{2};

if nargin > 3
    doplot = varargin{4};
else
    doplot = 0;
end

if nargin > 4
    distance = varargin{5};
else
    distance = 10; % points
end
if nargin > 2 % for detection of negative values
    arg = varargin{3}; % (arg = 2)  
else
    arg = 1;
end

inds = 0;
bin = gt(vec,thresh); % gives binary vector above thresh
diff_vec = squeeze(diff(bin));
inds = find(diff_vec > 0.8);
if arg == 2 % get also negative change
    inds = find(diff_vec ~= 0);
end
for i = 1:(length(inds)-1) % takes out close indices that are from the same stimulation
    if (inds(i+1) - inds(i)) < distance
        inds(i+1) =  inds(i);
    end
end
inds = unique(inds);
if length(inds) < 1
    inds = 0;
end

%% Checks ________________________________________

if doplot % add plots
    fg = 23; 
    figure(fg); clf;
    set(fg,'position',[ 1320         730         556         249]);
    set(fg,'name','findChange');
    subplot(2,1,1)
    % the stimulation trace
    plot(vec);
    hold on
    plot([0 length(vec)],[thresh thresh],'r')
    title('the input trace')
    
    subplot(2,1,2)
    % the diff and the treshold
    plot(diff_vec)
    hold on;
    if size(inds,2) > 1 % 
        inds = inds';
    end
    plot([inds  inds]',[zeros(1,length(inds)); ones(1,length(inds))],':g')
    title('The differentiated vector. The red line is the input theshold')
end
