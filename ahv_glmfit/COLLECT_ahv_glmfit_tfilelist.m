function [cfg_master, out, ts_pca1, ns_pca1] = COLLECT_ahv_glmfit_tfilelist(tfilelist, Z)
% JJS. 2024-10-02. This version works from a tfile list of .t files (path included), instead of pushing in and out of each session folder to calculate all neurons.
%                  This way we can calculate for only NPH neurons, for example.
% This version plots all points in the AHV and pupil tuning curves (in addition to the mean)
%% collect data from all sessions
%%
cfg_master = []; % overall params
cfg_master.doColor = 1;
cfg_master.FontSize = 10;
cfg_master.doPlot = 0;
cfg_master.WriteFig = 1;
cfg_master.dt = 0.005;
cfg_master.maxlag = 200; % bins for use in saccade PETH
cfg_master.debug = 0;
cfg_master.tc_binEdges = -150:10:150;
cfg_master.pupil_tc_binEdges = -80:5:80; % for pupilX TC

out = [];
for iNeuron = 1:length(tfilelist)
    [out{iNeuron}, ts_pca1(iNeuron), ns_pca1(iNeuron)] = SESSION_ahv_glmfit_tfilelist(cfg_master, tfilelist{iNeuron}, Z, iNeuron);
    
end

%% ERROR with concatenating. Look at matt's original version. 
% ahv_gain_fun = @(x) x.pca_sacc_both.rsq - x.pca_sacc.rsq;
% sacc_gain_fun = @(x) x.pca_sacc_both.rsq - x.ahv.rsq;
% rsq_fun = @(x) x.pca_sacc_both.rsq;
% 
% ahv_gain = arrayfun(ahv_gain_fun, out);
% sacc_gain = arrayfun(sacc_gain_fun, out);
% rsq_all = arrayfun(rsq_fun, out);

save('GLM.mat')

% PLOT
if cfg_master.doPlot == 1
    %% Overall Plot
    nCells = length(out);
    cfg_plot = [];
    cfg_plot.lim = 0.65;
    cfg_plot.rsq_bin = 0:0.05:1;
    cfg_plot.rsq_binC = cfg_plot.rsq_bin(1:end-1) + median(diff(cfg_plot.rsq_bin))/2;
    this_hist = histc(rsq_all, cfg_plot.rsq_bin);
    
    figure;
    
    subplot(222)           % bar plot of FULL model Rsq
    bar(cfg_plot.rsq_binC, this_hist(1:end-1));
    set(gca, 'TickDir', 'out', 'FontSize', 18);
    xlabel('full model R^2'); ylabel('count');
    
    subplot(221)           % Scatterplot of gain from AHV or EYE
    if cfg_master.doColor == 1; cmap = colormap(jet(round(100*max(rsq_all)))); else cmap = zeros(100,3); end
    for iC = 1:nCells
        h = text(ahv_gain(iC), sacc_gain(iC), num2str(iC));
        iColor = round(100*rsq_all(iC)); if iColor <= 0; iColor = 1; end
        set(h, 'FontSize', 16, 'FontWeight', 'Bold', 'Color', cmap(iColor,:));
        hold on;
    end
    plot([-0.1 cfg_plot.lim], [-0.1 cfg_plot.lim], 'k--');
    
    set(gca, 'XLim', [-0.1 cfg_plot.lim], 'YLim', [-0.1 cfg_plot.lim], 'TickDir', 'out', 'FontSize', 18);
    xlabel('R^2 gain from AHV'); ylabel('R^2 gain from eye movement');
    %% Inidvidual Tuning Curves
    cells_per_figure = 3; nFigures = ceil(nCells / cells_per_figure);
    xtight = .08;
    ytight = .04;
    savedir = ('D:\Jeff\U01\analysis\dot mat files\GLM\GLM subplots');
    
    for iF = 1:nFigures
        figure;
        start_cell = (iF-1)*cells_per_figure + 1;
        end_cell = min(start_cell + cells_per_figure-1, nCells);
        for iC = start_cell:end_cell
            plot_idx = iC - (iF-1)*cells_per_figure;
            %% subplot 1. AHV tuning curve.
            this_plot = (plot_idx-1)*3 + 1;
            %             subplot(cells_per_figure, 3, this_plot);
            subtightplot(cells_per_figure, 3, this_plot,[xtight ytight])
            hold on
            yyaxis left
            plot(out(iC).ahvscatterY', out(iC).ahvscatterX', '.')
            maxY = max(out(iC).ahvscatterX');
            yl = ylabel(num2str(iC)); set(yl, 'FontWeight', 'Bold');
            yyaxis right
            tc_bin_centers = cfg_master.tc_binEdges(1:end-1) + median(diff(cfg_master.tc_binEdges))/2;
            plot(tc_bin_centers, smoothdata(out(iC).tc), 'k', 'LineWidth', 5);
            box off;
            set(gca, 'TickDir', 'out', 'FontSize', cfg_master.FontSize, 'XLim', [-150 150], 'YLim', [0 maxY]);
            
            title(strcat('AHV only Rsq = '), num2str(out(iC).ahv.rsq)) %
            %%   subplot 2. Saccade PETH.
            %             subplot(cells_per_figure, 3, this_plot + 1);
            subtightplot(cells_per_figure, 3, this_plot + 1, [xtight ytight]);
            peth_bin_centers = -floor(length(out(iC).ns_peth)/2):floor(length(out(iC).ns_peth)/2);
            peth_bin_centers = peth_bin_centers .* cfg_master.dt;
            %peth_bin_centers = peth_bin_centers(1:end-1) + median(diff(peth_bin_centers))/2;
            plot(peth_bin_centers, out(iC).ns_peth ./ cfg_master.dt);
            hold on;
            plot(peth_bin_centers, out(iC).ts_peth ./ cfg_master.dt, 'r');
            set(gca, 'TickDir', 'out', 'FontSize', cfg_master.FontSize, 'XLim', [-.2 .2]); xlabel('');
            title(strcat('EYE only Rsq = '), num2str(out(iC).pca_sacc.rsq)) %
            %%  subplot 3. Pupil Tuning Curve.
            %             subplot(cells_per_figure, 3, this_plot + 2);
            subtightplot(cells_per_figure, 3, this_plot + 2, [xtight ytight]);
            tc_bin_centers = cfg_master.pupil_tc_binEdges(1:end-1) + median(diff(cfg_master.pupil_tc_binEdges))/2;
            plot(tc_bin_centers, out(iC).pupil_tc);
            set(gca, 'TickDir', 'out', 'FontSize', cfg_master.FontSize, 'XLim', [-60 60]); xlabel('');
            title(strcat('FULL Rsq = '), num2str(out(iC).pca_sacc_both.rsq)) %
        end
        if cfg_master.doPlot
            pushdir(savedir);
            WriteFig(num2str(iF))
            popdir;
        end
    end
end