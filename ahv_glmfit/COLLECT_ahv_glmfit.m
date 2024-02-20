%% collect data from all sessions
%% set up data path
cd('C:\data\U01\datatouse');
cfg = [];
%cfg.rats = {'M085', 'M089', 'M090'}; % only specific folders
f = dir; f = f(3:end); f = f([f.isdir]); % grab all folder names
cfg.rats = {f.name};

fd = getDataPath(cfg);

%%
cfg_master = []; % overall params
cfg_master.dt = 0.001;
cfg_master.maxlag = 200; % bins for use in saccade PETH
cfg_master.debug = 0;
cfg_master.tc_binEdges = -150:10:150;
cfg_master.pupil_tc_binEdges = -80:5:80; % for pupilX TC

out = [];

for iS = 1:length(fd)
   pushdir(fd{iS});
   
   %try
   out = cat(2, out, SESSION_ahv_glmfit(cfg_master));
   %catch
   %    disp('Session skipped.')
   %end
   
   popdir;
end

nCells = length(out);

%% plot
cfg_plot = [];
cfg_plot.lim = 0.5;
cfg_plot.rsq_bin = 0:0.05:1;

%ahv_gain_fun = @(x) x.both.rsq - x.sacc.rsq;
%sacc_gain_fun = @(x) x.both.rsq - x.ahv.rsq;
%rsq_fun = @(x) x.both.rsq;

ahv_gain_fun = @(x) x.pca_sacc_both.rsq - x.pca_sacc.rsq;
sacc_gain_fun = @(x) x.pca_sacc_both.rsq - x.ahv.rsq;
rsq_fun = @(x) x.pca_sacc_both.rsq;

ahv_gain = arrayfun(ahv_gain_fun, out);
sacc_gain = arrayfun(sacc_gain_fun, out);

figure;

subplot(221)
for iC = 1:nCells
    h = text(ahv_gain(iC), sacc_gain(iC), num2str(iC)); 
    set(h, 'FontSize', 16, 'FontWeight', 'Bold');
    hold on;
end
plot([-0.1 cfg_plot.lim], [-0.1 cfg_plot.lim], 'k--');

set(gca, 'XLim', [-0.1 cfg_plot.lim], 'YLim', [-0.1 cfg_plot.lim], 'TickDir', 'out', 'FontSize', 18);
xlabel('R^2 gain from AHV'); ylabel('R^2 gain from eye movement');

subplot(222)
rsq_all = arrayfun(rsq_fun, out);
this_hist = histc(rsq_all, cfg_plot.rsq_bin);

cfg_plot.rsq_binC = cfg_plot.rsq_bin(1:end-1) + median(diff(cfg_plot.rsq_bin))/2;
bar(cfg_plot.rsq_binC, this_hist(1:end-1));
set(gca, 'TickDir', 'out', 'FontSize', 18);
xlabel('full model R^2'); ylabel('count');

%%
cells_per_figure = 5; nFigures = ceil(nCells / cells_per_figure);

for iF = 1:nFigures
    
    figure;
    
    start_cell = (iF-1)*cells_per_figure + 1;
    end_cell = min(start_cell + cells_per_figure-1, nCells);
    
    for iC = start_cell:end_cell
        
        plot_idx = iC - (iF-1)*cells_per_figure;
        
        this_plot = (plot_idx-1)*3 + 1;
        subplot(cells_per_figure, 3, this_plot);
        
        tc_bin_centers = cfg_master.tc_binEdges(1:end-1) + median(diff(cfg_master.tc_binEdges))/2;
        plot(tc_bin_centers, out(iC).tc);
        box off;
        set(gca, 'TickDir', 'out', 'FontSize', 18, 'XLim', [-150 150]); xlabel('AHV');
        yl = ylabel(num2str(iC)); set(yl, 'FontWeight', 'Bold');
        
        subplot(cells_per_figure, 3, this_plot + 1);
        
        peth_bin_centers = -floor(length(out(iC).ns_peth)/2):floor(length(out(iC).ns_peth)/2);
        peth_bin_centers = peth_bin_centers .* cfg_master.dt;
        %peth_bin_centers = peth_bin_centers(1:end-1) + median(diff(peth_bin_centers))/2;
        
        plot(peth_bin_centers, out(iC).ns_peth ./ cfg_master.dt);
        hold on;
        plot(peth_bin_centers, out(iC).ts_peth ./ cfg_master.dt, 'r');
        set(gca, 'TickDir', 'out', 'FontSize', 18, 'XLim', [-1 1]); xlabel('time to sacc (s)');
        
        subplot(cells_per_figure, 3, this_plot + 2);
        
        tc_bin_centers = cfg_master.pupil_tc_binEdges(1:end-1) + median(diff(cfg_master.pupil_tc_binEdges))/2;
        plot(tc_bin_centers, out(iC).pupil_tc);
        set(gca, 'TickDir', 'out', 'FontSize', 18, 'XLim', [-60 60]); xlabel('pupil x position');
        
    end
end