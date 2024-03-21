%%
rng(0)
cd('C:\Jeff\U01\datatouse\M052\M052-2020-09-04');
sd = LoadSessionData([]);

cfg_master = [];
cfg_master.iCell = 1; % which cell to plot
cfg_master.nSaccades = 10; % if empty, plot all; if some number N, pick N random saccades

%%
cfg_peth = [];
cfg_peth.mode = 'interp';
cfg_peth.dt = 0.02;
cfg_peth.window = [-0.25 0.25];
[tsd] = TSDpeth(cfg_peth, sd.tsdH, sd.temporalSaccades);
% [p, a] = TSDpeth(cfg_peth, sd.tsdH, sd.temporalSaccades);
% a = tsd; 

% select which saccades to keep
N = size(a.data, 1);
if ~isempty(cfg_master.nSaccades)
    keep_idx = randperm(N);
    keep_idx = keep_idx(1:cfg_master.nSaccades);
else
    keep_idx = 1:N;
end

subplot(221);
plot(tsd.tvec, tsd.data(keep_idx,:), 'Color', [0.7 0.7 0.7]);
hold on;
plot(p.tvec, p.data, 'k', 'LineWidth', 2)
xline(0, ':', 'Color', 'r', 'LineWidth', 2);

set(gca, 'XLim', cfg_peth.window, 'TickDir', 'out', 'FontSize', 24, 'XTick', [-0.25:0.125:0.25], 'XTickLabel', {'-0.25', '', '0', '', '0.25'}); box off;
ylabel('horiz. eye pos (a.u.)');


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

for iSess = 1:length(fdr);
    [a{iSess} b{iSess} c{iSess}] = fileparts(fdr); 
end
%% plot the basic data 

numtemporalSaccadesManual
numtemporalSaccadesNew = length(temporalSaccadesNew)
temporal_overlap = intersect(temporalSaccadesManual, temporalSaccadesNew);
num_temporal_overlap = length(temporal_overlap)
num_temporal_miss = abs(length(temporalSaccadesManual) - length(temporalSaccadesNew))

numnasalSaccadesManual
numnasalSaccadesNew = length(nasalSaccadesNew)
nasal_overlap = intersect(nasalSaccadesManual, nasalSaccadesNew);
num_nasal_overlap = length(nasal_overlap)
num_nasal_miss = abs(length(nasalSaccadesManual) - length(nasalSaccadesNew))


clf
hold on
SSN = HD_GetSSN;

        plot(tsdH.tvec, tsdH.data, 'Color', [.301 .745 .933], 'LineStyle', '--') % HORIZONTAL PUPIL POSITION 
        plot(diffH.tvec, diffH.data, 'Color', [.85 .325 .098], 'LineStyle', '-') % HORIZONTAL PUPIL VELOCITY 
        plot(diffV.tvec, diffV.data, 'm', 'LineStyle', '-') % VERTICAL PUPIL VELOCITY 


        plot(temporalSaccadesManual, temporalAmplitudesManual, 'r.', 'MarkerSize', 25)
        plot(nasalSaccadesManual, nasalAmplitudesManual, 'r.', 'MarkerSize', 25)
        
        plot(temporalSaccadesNew, temporalAmplitudesNew, 'go', 'MarkerSize', 25)
        plot(nasalSaccadesNew, nasalAmplitudesNew, 'go', 'MarkerSize', 25)

%         line([c(1) c(2)], [-3 -3], 'Color', 'r', 'LineWidth', 3, 'LineStyle', '--', 'Color', 'g')
%         line([c(1) c(2)], [-3 -3], 'Color', 'r', 'LineWidth', 3, 'LineStyle', '--', 'Color', 'g')
        c = axis;

        line([c(1) c(2)], [10.02 10.02], 'Color', 'k', 'LineWidth', 3, 'LineStyle', '--')
        line([c(1) c(2)], [-10.02 -10.02], 'Color', 'k', 'LineWidth', 3, 'LineStyle', '--')

        line([c(1) c(2)], [cfg_in.threshT cfg_in.threshT], 'Color', 'g', 'LineWidth', 3, 'LineStyle', '--')
        line([c(1) c(2)], [cfg_in.threshN cfg_in.threshN], 'Color', 'g', 'LineWidth', 3, 'LineStyle', '--')
        
        yyaxis right
        plot(AHV_tsd.tvec, AHV_tsd.data, 'Color', [.75 .75 0])
        
        set(gca, 'FontSize', 20)
        title(SSN)




    



















