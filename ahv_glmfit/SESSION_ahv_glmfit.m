function out = SESSION_ahv_glmfit(cfg_in)
% function out = SESSION_ahv_glmfit(cfg_in)
%
% perform single session GLM fits of spiking data with AHV and saccade
% regressors
%
% run from single session's data folder
%

cfg_def = []; % overall params
cfg_def.dt = 0.001;
cfg_def.maxlag = 500; % bins for use in saccade PETH
cfg_def.debug = 0;
cfg_def.tc_binEdges = -150:10:150;

cfg = ProcessConfig2(cfg_def, cfg_in);

cfg_master = cfg;
cfg_master.gausswin = gausswin(1./cfg_master.dt, 20); cfg_master.gausswin = cfg_master.gausswin ./ sum(cfg_master.gausswin);

%% load data
AHV_tsd = Get_AHV([]);
S = LoadSpikesJeff; nCells = length(S.t);
[temporalSaccades, nasalSaccades, ~, ~, ~, pupilX_tsd, pupilY_tsd] = processPupilData2([]);

%%
p = table; % contains regressors

nMaxVars = 5; % only used for initializing t-stat matrix
% baseline model MUST be defined first or things will break!
sd.m.baseline.modelspec = 'spk ~ 1 + time';
%sd.m.ahv.modelspec = 'spk ~ 1 + time + ns';
sd.m.ahv.modelspec = 'spk ~ 1 + time + ahv';
sd.m.sacc.modelspec = 'spk ~ 1 + time + ns + ts';
sd.m.both.modelspec = 'spk ~ 1 + time + ns + ts + ahv';

% init error vars
mn = fieldnames(sd.m);
for iM = 1:length(mn)
   %sd.m.(mn{iM}).err = zeros(nCells, length(sd.TVECc)); % needs to be zeros because error output will be added to this
   sd.m.(mn{iM}).tstat = nan(nCells, nMaxVars);
end

% common timebase for this session
TVECe = 0:cfg_master.dt:AHV_tsd.tvec(end); % edges
sd.TVECc = TVECe(1:end-1)+cfg_master.dt/2; % centers

p.time = sd.TVECc'; % elapsed time predictor

p.ahv = interp1(AHV_tsd.tvec, AHV_tsd.data, sd.TVECc)'; % AHV predictor

% saccade predictors
[ns_binarized, ns_bin] = histc(nasalSaccades, TVECe); ns_binarized = ns_binarized(1:end - 1);
[ts_binarized, ts_bin] = histc(temporalSaccades, TVECe); ts_binarized = ts_binarized(1:end - 1);

% PCA PETH predictors
% pushdir('C:\Users\mvdm\Dropbox\projects\Jeff\ahv_peth_pca');
% load SaccadePETHs
% popdir;
% 
% x1 = FRxBinNsmooth'; x2 = FRxBinTsmooth'; x = cat(2, x1, x2);
% [coeff, score] = pca(x);
% 
% pca1_kernel = interp1(binCenters, score(:, 1), binCenters(1):cfg_master.dt:binCenters(end));
% pca1_kernel = pca1_kernel .* gausswin(length(pca1_kernel), 3)';
% p.ts_pca1 = conv2(ts_binarized, pca1_kernel, 'same')';
% p.ns_pca1 = conv2(ns_binarized, pca1_kernel, 'same')';
% 
% pca2_kernel = interp1(binCenters, score(:, 2), binCenters(1):cfg_master.dt:binCenters(end));
% pca2_kernel = pca2_kernel .* gausswin(length(pca2_kernel), 3)';
% p.ts_pca2 = conv2(ts_binarized, pca2_kernel, 'same')';
% p.ns_pca2 = conv2(ns_binarized, pca2_kernel, 'same')';
% 
% pca3_kernel = interp1(binCenters, score(:, 3), binCenters(1):cfg_master.dt:binCenters(end));
% pca3_kernel = pca3_kernel .* gausswin(length(pca3_kernel), 3)';
% p.ts_pca3 = conv2(ts_binarized, pca3_kernel, 'same')';
% p.ns_pca3 = conv2(ns_binarized, pca3_kernel, 'same')';

%% loop over cells
cellCount = 0; 
for iC = nCells:-1:1

    fprintf('Cell %d/%d...\n',iC,nCells);
    
    % figure out if we should include this cell
    
    % dependent variable: binned spike train
    spk_binned = histc(S.t{iC}, TVECe); spk_binned = spk_binned(1:end - 1);
    spk_binned = logical(spk_binned > 0); % binarize
    
    p.spk = spk_binned;
    p.spk = conv2(p.spk, cfg_master.gausswin, 'same') ./ cfg_master.dt;
    
    % get AHV tc
    this_ahv = tsd(sd.TVECc, p.ahv'); % could be done outside loop
    this_sdf = tsd(sd.TVECc, p.spk');
    
    cfg_tc = [];
    cfg_tc.binEdges = cfg_master.tc_binEdges;
    cfg_tc.dt = median(diff(sd.TVECc));
    tc = TuningCurvesSDF(cfg_tc, this_sdf, this_ahv);
    out(iC).tc = tc.tc;
    
    % saccade predictors based on PETHs
    this_xc = xcorr(p.spk, ns_binarized, cfg_master.maxlag*2, 'unbiased');
    out(iC).ns_peth = this_xc;
        
    this_xc = this_xc .* gausswin(length(this_xc), 3); % smooth edges, works better as regressor
    ns_conv = conv2(ns_binarized, this_xc', 'same');
    %ns_conv = conv2(ns_binarized, this_xc'-mean(this_xc), 'same');
    %ns_conv = conv2(ns_conv, cfg_master.gausswin', 'same');
    
    this_xc = xcorr(p.spk, ts_binarized, cfg_master.maxlag*2, 'unbiased');
    out(iC).ts_peth = this_xc;
    
    this_xc = this_xc .* gausswin(length(this_xc), 3);
    ts_conv = conv2(ts_binarized, this_xc', 'same');
    %ts_conv = conv2(ts_binarized, this_xc'-mean(this_xc), 'same');
    %ts_conv = conv2(ts_conv, cfg_master.gausswin', 'same');
   
    p.ns = ns_conv';
    p.ts = ts_conv';
    p.sacc = ns_conv' + ts_conv';
    
    %this_m = fitglm(p, sd.m.(mn{4}).modelspec, 'Distribution', 'binomial')
    for iM = 1:length(mn)
        fprintf('Model %s:\n', mn{iM});
        this_m = fitglm(p, sd.m.(mn{iM}).modelspec, 'Distribution', 'poisson')
        
        out(iC).(mn{iM}).rsq = this_m.Rsquared.Adjusted;
    end
    
    if cfg.debug
        this_fit = this_m.predict(p);
           
        plot(sd.TVECc, p.spk);
        hold on;
        plot(sd.TVECc, this_fit, 'r');
        plot(sd.TVECc(ns_bin), 0, '.c', 'MarkerSize', 20);
        plot(sd.TVECc(ts_bin), 0, '.m', 'MarkerSize', 20);
        set(gca, 'FontSize', 18); ylabel('Firing rate (Hz)'); xlabel('time (s)');
        
        pause;
    end
    
end