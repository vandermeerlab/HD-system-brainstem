%%
rng(0)
% cd('C:\Jeff\U01\datatouse\M052\M052-2020-09-04');
sd = LoadSessionData([]);

cfg_master = [];
cfg_master.iCell = 1; % which cell to plot
cfg_master.nSaccades = 10; % if empty, plot all; if some number N, pick N random saccades

%% Top Left Subplot. Do TEMPORAL saccades first 
cfg_peth = [];
cfg_peth.mode = 'interp';
cfg_peth.dt = 0.02;
cfg_peth.window = [-0.25 0.25];
[p, a] = TSDpeth(cfg_peth, sd.tsdH, sd.temporalSaccades);

% select which saccades to keep
N = size(a.data, 1);
if ~isempty(cfg_master.nSaccades)
    keep_idx = randperm(N);
    keep_idx = keep_idx(1:cfg_master.nSaccades);
else
    keep_idx = 1:N;
end

subplot(221);
plot(a.tvec, a.data(keep_idx,:), 'Color', [0.7 0.7 0.7]);
hold on;
plot(p.tvec, p.data, 'k', 'LineWidth', 2)
xline(0, ':', 'Color', 'r', 'LineWidth', 2);

set(gca, 'XLim', cfg_peth.window, 'TickDir', 'out', 'FontSize', 24, 'XTick', [-0.25:0.125:0.25], 'XTickLabel', {'-0.25', '', '0', '', '0.25'}); box off;
ylabel('horiz. eye pos (a.u.)');
% title(sd.fn(cfg_master.iCell))
[a, b, c,] = fileparts(sd.fc{cfg_master.iCell});
cellID = strcat(b,c);
title(cellID);
%%
this_S = SelectTS([], sd.S, cfg_master.iCell);

cfg_peth.doPlot = 0;
[outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm(cfg_peth, this_S, sd.temporalSaccades);

subplot(223)
for iT = 1:length(keep_idx)

    this_idx = find(outputT == keep_idx(iT)); % spikePETH idxs for this trial
    plot(outputS(this_idx), iT, '|k', 'MarkerSize', 20); hold on;

end
set(gca, 'XLim', cfg_peth.window, 'TickDir', 'out', 'FontSize', 24, 'XTick', [-0.25:0.125:0.25], 'XTickLabel', {'-0.25', '', '0', '', '0.25'}, 'YTick', [], 'YColor', [1 1 1]); box off;
ax = gca; ax.YAxis.Visible = 'off'; ax.Color = 'none'
xlabel('time (s)');


if strcmp(sd.ExpKeys.TTHemisphere{cfg_master.iCell}, {'L'})
    title('IPSILATERAL SACCADES') 
elseif strcmp(sd.ExpKeys.TTHemishphere{cfg_master.iCell}, {'R'})
    title('CONTRALATERAL SACCADES') 
else 
    warning('hemisphere unknown')
end

%%
[p, a] = TSDpeth(cfg_peth, sd.tsdH, sd.nasalSaccades);

% select which saccades to keep
N = size(a.data, 1);
if ~isempty(cfg_master.nSaccades)
    keep_idx = randperm(N);
    keep_idx = keep_idx(1:cfg_master.nSaccades);
else
    keep_idx = 1:N;
end

subplot(222);
plot(a.tvec, a.data(keep_idx,:), 'Color', [0.7 0.7 0.7]);
hold on;
plot(p.tvec, p.data, 'k', 'LineWidth', 2)
xline(0, ':', 'Color', 'r', 'LineWidth', 2);

set(gca, 'XLim', cfg_peth.window, 'TickDir', 'out', 'FontSize', 24, 'XTick', [-0.25:0.125:0.25], 'XTickLabel', {'-0.25', '', '0', '', '0.25'}); box off;
ylabel('horiz. eye pos (a.u.)');

% title(sd.fn(cfg_master.iCell))
[a, b, c,] = fileparts(sd.fc{cfg_master.iCell});
cellID = strcat(b,c);
title(cellID);

%%
this_S = SelectTS([], sd.S, cfg_master.iCell);

cfg_peth.doPlot = 0;    
[outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm(cfg_peth, this_S, sd.nasalSaccades);

subplot(224)
for iT = 1:length(keep_idx)

    this_idx = find(outputT == keep_idx(iT)); % spikePETH idxs for this trial
    plot(outputS(this_idx), iT, '|k', 'MarkerSize', 20); hold on;

end
set(gca, 'XLim', cfg_peth.window, 'TickDir', 'out', 'FontSize', 24, 'XTick', [-0.25:0.125:0.25], 'XTickLabel', {'-0.25', '', '0', '', '0.25'}, 'YTick', [], 'YColor', [1 1 1]); box off;
ax = gca; ax.YAxis.Visible = 'off'; ax.Color = 'none'
xlabel('time (s)');

if strcmp(sd.ExpKeys.TTHemisphere{cfg_master.iCell}, {'R'})
    title('IPSILATERAL SACCADES') 
elseif strcmp(sd.ExpKeys.TTHemisphere{cfg_master.iCell}, {'L'})
    title('CONTRALATERAL SACCADES') 
else 
    warning('hemisphere unknown')
end

