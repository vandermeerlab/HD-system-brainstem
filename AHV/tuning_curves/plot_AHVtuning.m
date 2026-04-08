%% cd to data folder
% miscellaney
FontSize = 16;
histXmin = 0.01;  
histXmax = 0.2;
SSN = HD_GetSSN('SingleSession');
%% get AHV
cfg_AHV = [];
cfg_AHV.subsample_factor = 10;
[AHV_tsd, tc_out] = AHV_tuning(cfg_AHV, S);
% [AHV_tsd, S, tc_out] = AHV_tuning(cfg_AHV);
AHV_dt = median(diff(AHV_tsd.tvec));

%% plot tuning curves
fc = FindFiles('*.t', 'CheckSubdirs', 0);
numCells = size(tc_out.tc, 1);
clf

for iCell = 1:numCells; % first row is Tuning Curves
    [a, b, c] = fileparts(fc{iCell});
    subplot(5,numCells,iCell)
    set(gca, 'TickDir', 'out', 'FontSize', FontSize)
    plot(tc_out.usr.binCenters, tc_out.tc(iCell,:), 'MarkerSize', 5);
    xlabel('AHV (deg./sec)', 'FontSize', FontSize)
    ylabel('FR (Hz)', 'FontSize', FontSize)
    set(groot, 'DefaultLegendInterpreter', 'none')
    title(strcat(SSN, '-', b), 'FontSize', FontSize)
    
end

%% plot scatterplot
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = AHV_dt;
F = MakeQfromS(cfg_Q, S);

% convert to FR
F.data = F.data ./ cfg_Q.dt;

% find FR corresponding to each AHV sample
F_idx = nearest_idx3(AHV_tsd.tvec, F.tvec);
AHV_F = F.data(:,F_idx);

iCell = 0;
for iPlot = numCells+1: 2*numCells;   % second row is Firing Rate scatterplot
    iCell = iCell +1;
    subplot(5,numCells,iPlot)
    set(gca, 'TickDir', 'out', 'FontSize', FontSize)
    plot(AHV_tsd.data, AHV_F(iCell,:), '.', 'MarkerSize', .5);
    xlabel('AHV (deg./sec)', 'FontSize', FontSize)
    ylabel('FR (Hz)', 'FontSize', FontSize)
end

%% acf
cfg_acf = [];
cfg_acf.binsize = 0.001;
cfg_acf.max_t = 0.5;
cfg_acf.smooth = 1; % set to 1 to compute ccf on SDF, 0 for raw spikes
cfg_acf.gauss_w = 1; % width of Gaussian convolution window (in s)
cfg_acf.gauss_sd = 0.005; % SD of Gaussian convolution window (in s)

iCell = 0;
for iPlot = 2*numCells+1: 3*numCells;   % third row is the autocorrelation
    iCell = iCell +1;
    subplot(5,numCells,iPlot)
    set(gca, 'TickDir', 'out', 'FontSize', FontSize)
    [acf, tvec] = ccf(cfg_acf, S.t{iCell}, S.t{iCell});
    midpoint = ceil(length(acf)./2);
    acf(midpoint) = 0;
    plot(tvec, acf);
    xlabel('Time (sec)', 'FontSize', FontSize)
    set(gca, 'XTick', -0.5:0.1:0.5); grid on;
end

%% HistISI
iCell = 0;
for iPlot = 3*numCells+1: 4*numCells;   % fourth row is the HistISI
    iCell = iCell +1;
    subplot(5,numCells,iPlot); hold on
%     set(gca, 'TickDir', 'out', 'FontSize', FontSize)
    [h, n] = HistISIsubplot(S.t{iCell});
    HistISIsubplot(S.t{iCell});
    [c, i] = max(h);
    line([n(i) n(i)], [0 162]);
    
    grid on 
    set(gca, 'TickDir', 'out', 'XLim', [histXmin histXmax], 'FontSize', FontSize)
    xlabel('Time (sec)', 'FontSize', FontSize)
end

%% Firing Rate  
iCell = 0;
for iPlot = 4*numCells+1: 5*numCells;   % firing rate plot for the whole session
    iCell = iCell +1;
    subplot(5,numCells,iPlot)
    set(gca, 'TickDir', 'out', 'XLim', [0.01 0.5])
    cfg_Q = []; cfg_Q.dt = 0.001; cfg_Q.gausswin_sd = 0.05;cfg_Q.smooth = 'gauss';
    Q = MakeQfromS(cfg_Q, S);
    plot(Q.tvec, Q.data./cfg_Q.dt)
    
    
    
end




