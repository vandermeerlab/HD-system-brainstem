%% collector script to run GLM fits on all sessions, collect and plot the results
% set up list of sessions: run from top-level data folder, e.g. cd('C:\data\U01\datatouse');
cfg = [];
cfg.rats = {'M039', 'M052', 'M055', 'M079', 'M080', 'M085', 'M086', 'M089', 'M090', 'M094', 'M096', 'M104', 'M105', 'M112', 'M212', 'M269', 'M271', 'M293'};
%cfg.rats = {'M052'};
fd = getDataPath(cfg);

endSess = length(fd);
%%
cfg_master = []; % overall params
cfg_master.dt = 0.001;
cfg_master.maxlag = 200; % bins for use in saccade PETH
cfg_master.debug = 0;
cfg_master.tc_binEdges = -150:10:150;
cfg_master.pupil_tc_binEdges = -80:5:80; % for pupilX TC
cfg_master.doPlot = 1;
cfg_master.doColor = 1; 
cfg_master.FontSize = 10;
cfg_master.intermediate_file_path = 'C:\data\U01\datatouse\intermediate_files';

out = [];
sessUsed = {};
sessCounter = 0;
skip = 0;  % there are four sessions that return errors when the model that includes saccade amplitude is used.

%process_varargin(varargin);

for iS = 1:endSess
    pushdir(fd{iS});
    SSN = HD_GetSSN; disp(SSN);
    if exist(strcat(SSN, '-VT1_proc.mat')) == 2
        if skip == 1
            if strcmp(SSN, 'M094-2021-12-28') == 0 && strcmp(SSN, 'M105-2021-01-17') == 0 && strcmp(SSN, 'M212-2021-07-21') == 0 && strcmp(SSN, 'M271-2021-08-28') == 0 && strcmp(SSN, 'M086-2020-11-21') == 0 && strcmp(SSN, 'M094-2020-12-28') == 0
                sessCounter = sessCounter +1;
                sessUsed{sessCounter} = SSN;
                out = cat(2, out, SESSION_ahv_glmfit_jeff(cfg_master));
            else
                disp('problem session. skipping for now...')
            end
        else
            sessCounter = sessCounter +1;
            sessUsed{sessCounter} = SSN;
            out = cat(2, out, SESSION_ahv_glmfit_jeff(cfg_master));
        end
    else
        disp('no video tracking data. skipping session...')
    end
    popdir;
end
%% PLOT
if cfg_master.doPlot == 1
    nCells = length(out);
    %% plot
    cfg_plot = [];
    cfg_plot.lim = 0.65;
    cfg_plot.rsq_bin = 0:0.05:1;
    
    ahv_gain_fun = @(x) x.pca_sacc_both.rsq - x.pca_sacc.rsq;
    sacc_gain_fun = @(x) x.pca_sacc_both.rsq - x.ahv.rsq;
    rsq_fun = @(x) x.pca_sacc_both.rsq;
    
    ahv_gain = arrayfun(ahv_gain_fun, out);
    sacc_gain = arrayfun(sacc_gain_fun, out);
    
    figure;
    
    subplot(222)           % bar plot of FULL model Rsq
    rsq_all = arrayfun(rsq_fun, out);
    this_hist = histc(rsq_all, cfg_plot.rsq_bin);
    
    cfg_plot.rsq_binC = cfg_plot.rsq_bin(1:end-1) + median(diff(cfg_plot.rsq_bin))/2;
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
    %%
    cells_per_figure = 5; nFigures = ceil(nCells / cells_per_figure);
    for iF = 1:nFigures
        figure;
        start_cell = (iF-1)*cells_per_figure + 1;
        end_cell = min(start_cell + cells_per_figure-1, nCells);
        for iC = start_cell:end_cell
            plot_idx = iC - (iF-1)*cells_per_figure;
            %% subplot 1. AHV tuning curve.
            this_plot = (plot_idx-1)*3 + 1;
            subplot(cells_per_figure, 3, this_plot);
            tc_bin_centers = cfg_master.tc_binEdges(1:end-1) + median(diff(cfg_master.tc_binEdges))/2;
            plot(tc_bin_centers, out(iC).tc);
            %             subtightplot(cells_per_figure, 3, this_plot,[tc_bin_centers, out(iC).tc])
            box off;
            set(gca, 'TickDir', 'out', 'FontSize', cfg_master.FontSize, 'XLim', [-150 150]); xlabel('AHV');
            yl = ylabel(num2str(iC)); set(yl, 'FontWeight', 'Bold');
            title(strcat('AHV only Rsq = '), num2str(out(iC).ahv.rsq)) %
            %%   subplot 2. Saccade PETH.
            subplot(cells_per_figure, 3, this_plot + 1);
            peth_bin_centers = -floor(length(out(iC).ns_peth)/2):floor(length(out(iC).ns_peth)/2);
            peth_bin_centers = peth_bin_centers .* cfg_master.dt;
            %peth_bin_centers = peth_bin_centers(1:end-1) + median(diff(peth_bin_centers))/2;
            plot(peth_bin_centers, out(iC).ns_peth ./ cfg_master.dt);
            hold on;
            plot(peth_bin_centers, out(iC).ts_peth ./ cfg_master.dt, 'r');
            set(gca, 'TickDir', 'out', 'FontSize', cfg_master.FontSize, 'XLim', [-.2 .2]); xlabel('time to sacc (s)');
            title(strcat('EYE only Rsq = '), num2str(out(iC).pca_sacc.rsq)) %
            %%  subplot 3. Pupil Tuning Curve.
            subplot(cells_per_figure, 3, this_plot + 2);
            tc_bin_centers = cfg_master.pupil_tc_binEdges(1:end-1) + median(diff(cfg_master.pupil_tc_binEdges))/2;
            plot(tc_bin_centers, out(iC).pupil_tc);
            set(gca, 'TickDir', 'out', 'FontSize', cfg_master.FontSize, 'XLim', [-60 60]); xlabel('pupil x position');
            title(strcat('FULL Rsq = '), num2str(out(iC).pca_sacc_both.rsq)) %
        end
    end
end