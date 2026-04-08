function [out, ahv_gain, sacc_gain, rsq_all, ts_pca1, ns_pca1] = COLLECT_ahv_glmfit(startSess, endSess, fd)
% original function by MvdM.
% 2024-10-22. ***Edited by JJS to skip sessions where there is an error [~strcmp(pwd, __)]. These sessions won't run b/c in SESSION_ahv_glmfit.m
% the size of ns_idx or ts_idx is n-1 the size of nasalAmplitudes or temporalAmplitudes, respectively. Need to figure out why this is happening
% and fix. I don't remember this cropping up before. This function works through the entire directory. Future version will take a list of tfiles.

% Leave startSess and endSess as empty ([]) if you want to start at folder 1 and proceed through to the end.
%% collect data from all sessions
%% set up data path
% cd('C:\Jeff\U01\datatouse');
% cfg = [];
% %cfg.rats = {'M085', 'M089', 'M090'}; % only specific folders
% f = dir; f = f(3:end); f = f([f.isdir]); % grab all folder names
% cfg.rats = {f.name};
% fd = getDataPath(cfg);

% if isempty(fd)
%     fd = FindFiles('*keys.m');
% end

%%
cfg_master = []; % overall params
cfg_master.dt = 0.001;
cfg_master.maxlag = 200; % bins for use in saccade PETH
cfg_master.debug = 0;
cfg_master.tc_binEdges = -150:10:150;
cfg_master.pupil_tc_binEdges = -80:5:80; % for pupilX TC

out = [];
ts_pca1 = [];
ns_pca1 = [];

if isempty(startSess); startSess = 1; end
if isempty(endSess); endSess = length(fd); end
if isempty(fd); error('not sessions given as input'); end

for iS = startSess : endSess
    pushdir(fd{iS});
    if ~strcmp(pwd, 'C:\Jeff\U01\datatouse\M086\M086-2020-11-21') && ~strcmp(pwd, 'C:\Jeff\U01\datatouse\M094\M094-2020-12-28') ...
            && ~strcmp(pwd, 'C:\Jeff\U01\datatouse\M105\M105-2021-01-17') && ~strcmp(pwd, 'C:\Jeff\U01\datatouse\M212\M212-2021-07-21') ...
            && ~strcmp(pwd, 'C:\Jeff\U01\datatouse\M222\M222-2022-01-23') && ~strcmp(pwd, 'C:\Jeff\U01\datatouse\M271\M271-2021-08-28') ...
            && ~strcmp(pwd, 'C:\Jeff\U01\datatouse\M279\M279-2022-06-21-2') && ~strcmp(pwd, 'C:\Jeff\U01\datatouse\M282\M282-2022-02-04-1') ...
            && ~strcmp(pwd, 'C:\Jeff\U01\datatouse\M287\M287-2022-08-05') && ~strcmp(pwd, 'C:\Jeff\U01\datatouse\M287\M287-2022-08-11-1')   ...
            && ~strcmp(pwd, 'C:\Jeff\U01\datatouse\M293\M293-2021-08-19') && ~strcmp(pwd, 'C:\Jeff\U01\datatouse\M391\M391-2022-12-12') ...
            && ~strcmp(pwd, 'C:\Jeff\U01\datatouse\M402\M402-2022-11-02-2')
        % error in these sessions
        %try
        %         out = cat(2, out, SESSION_ahv_glmfit(cfg_master));
        [this_out, this_ts_pca1, this_ns_pca1] = SESSION_ahv_glmfit(cfg_master);
        
        out = cat(2, out, this_out);
        
        ts_pca1 = cat(2, this_ts_pca1, ts_pca1);
        ns_pca1 = cat(2, this_ns_pca1, ns_pca1);
        %catch
        %    disp('Session skipped.')
        %end
    else
        disp('skipping session.................................................................................................')
    end
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