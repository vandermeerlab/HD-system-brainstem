function [out, ts_pca1, ns_pca1] = SESSION_ahv_glmfit_tfilelist(cfg_in, tfile, Z, iNeuron)
% function out = SESSION_ahv_glmfit(cfg_in)
%
% perform single session GLM fits of spiking data with AHV and saccade
% regressors
%
% run from single session's data folder
%
%      Inputs:
%                   Z - structure with saccade peth matrix for all NPH neurons

cfg_def = []; % overall params
cfg_def.dt = 0.005;
cfg_def.maxlag = 500; % bins for use in saccade PETH
cfg_def.debug = 1;
cfg_def.tc_binEdges = -150:10:150; % for AHV TC
cfg_def.pupil_tc_binEdges = -80:5:80; % for pupilX TC
cfg_def.skip_amp = 1;  % skip calculating saccade amplitude as part of the model
cfg_def.plotAllPoints = 1;   % plot all of the points in the tuning curve scatterplot and all of the points in the pupil position plot
cfg_def.intermediate_file_path = 'C:\Jeff\U01\datatouse'; % needed for SaccadePETH.mat file, define this in function call

cfg_master = ProcessConfig2(cfg_def, cfg_in);
cfg_master.gausswin = gausswin(1./cfg_master.dt, 20); cfg_master.gausswin = cfg_master.gausswin ./ sum(cfg_master.gausswin);

[path, neuron_to_use, ext] = fileparts(tfile);
neuronID = strcat(neuron_to_use, ext);
disp(neuron_to_use)

pushdir(path);
SSN = HD_GetSSN;
EvalKeys

%% load data
AHV_tsd = Get_AHV([]);
cfg_spikes.fc = {neuronID};
S = LoadSpikesJeff(cfg_spikes); 
assert(length(S.t) == 1);

load(FindFile('*saccades-edited.mat'))
keep = find(~isnan(temporalSaccades));
temporalSaccades = temporalSaccades(keep); temporalAmplitudes = temporalAmplitudes(keep);

keep = find(~isnan(nasalSaccades));
nasalSaccades = nasalSaccades(keep); nasalAmplitudes = nasalAmplitudes(keep);

%%
p = table; % contains regressors

nMaxVars = 5; % only used for initializing t-stat matrix
% % baseline model MUST be defined first or things will break!
% sd.m.baseline.modelspec = 'spk ~ 1 + time';
% %sd.m.ahv.modelspec = 'spk ~ 1 + time + ns';
% % sd.m.ahv.modelspec = 'spk ~ 1 + time + ahv';
% %sd.m.sacc.modelspec = 'spk ~ 1 + time + ns + ts';
% %sd.m.both.modelspec = 'spk ~ 1 + time + ns + ts + ahv';
% sd.m.pca_sacc.modelspec = 'spk ~ 1 + time + ts_pca1 + ns_pca1 + ts_pca2 + ns_pca2 + ts_pca3 + ns_pca3 + pupilX + pupilV';
% sd.m.pca_sacc_both.modelspec = 'spk ~ 1 + time + ts_pca1 + ns_pca1 + ts_pca2 + ns_pca2 + ts_pca3 + ns_pca3 + pupilX + pupilV + ahv';
% %sd.m.pca_sacc.modelspec = 'spk ~ 1 + time + pupilX + ts_pca1*ts_ampl + ns_pca1*ns_ampl + ts_pca2*ts_ampl + ns_pca2*ns_ampl + ts_pca3*ts_ampl + ns_pca3*ns_ampl';
% %sd.m.pca_sacc_both.modelspec = 'spk ~ 1 + time + pupilX + ahv + ts_pca1*ts_ampl + ns_pca1*ns_ampl + ts_pca2*ts_ampl + ns_pca2*ns_ampl + ts_pca3*ts_ampl + ns_pca3*ns_ampl';


sd.m.baseline.modelspec = 'spk ~ 1 + time';
sd.m.pca_sacc.modelspec = 'spk ~ 1 + time + ts_pca1 + ns_pca1 + ts_pca2 + ns_pca2 + ts_pca3 + ns_pca3';
% sd.m.pca_sacc_both.modelspec = 'spk ~ 1 + time + ts_pca1 + ns_pca1 + ts_pca2 + ns_pca2 + ts_pca3 + ns_pca3 + pupilX + pupilV + ahv';


% init error vars
mn = fieldnames(sd.m);
for iM = 1:length(mn)
    %sd.m.(mn{iM}).err = zeros(nCells, length(sd.TVECc)); % needs to be zeros because error output will be added to this
    sd.m.(mn{iM}).tstat = nan(1, nMaxVars);
end

% common timebase for this session
TVECe = 0:cfg_master.dt:AHV_tsd.tvec(end); % edges
sd.TVECc = TVECe(1:end-1)+cfg_master.dt/2; % centers

p.time = sd.TVECc'; % elapsed time predictor

p.ahv = interp1(AHV_tsd.tvec, AHV_tsd.data, sd.TVECc)'; % AHV predictor
this_ahv = tsd(sd.TVECc, p.ahv'); % used for making TCs within cell loop later

p.pupilX = interp1(tsdH.tvec, tsdH.data, sd.TVECc)'; % pupil X-position predictor
this_pupilX = tsd(sd.TVECc, p.pupilX');

p.pupilV = interp1(tsdH.tvec(2:end), diff(tsdH.data), sd.TVECc)';
this_pupilV = tsd(sd.TVECc, p.pupilV');

% saccade predictors
[ns_binarized, ns_bin] = histc(nasalSaccades, TVECe); ns_binarized = ns_binarized(1:end - 1);
[ts_binarized, ts_bin] = histc(temporalSaccades, TVECe); ts_binarized = ts_binarized(1:end - 1);

% PCA PETH predictors
% pushdir(cfg_master.intermediate_file_path); load SaccadePETHs; popdir;

x1 = Z.FRxBinNsmooth'; x2 = Z.FRxBinTsmooth'; x = cat(2, x1, x2);
[coeff, score] = pca(x);

pca1_kernel = interp1(Z.binCenters, score(:, 1), Z.binCenters(1): cfg_master.dt: Z.binCenters(end));
pca1_kernel = pca1_kernel .* gausswin(length(pca1_kernel), 3)';
p.ts_pca1 = conv2(ts_binarized, pca1_kernel, 'same')';
p.ns_pca1 = conv2(ns_binarized, pca1_kernel, 'same')';

pca2_kernel = interp1(Z.binCenters, score(:, 2), Z.binCenters(1): cfg_master.dt: Z.binCenters(end));
pca2_kernel = pca2_kernel .* gausswin(length(pca2_kernel), 3)';
p.ts_pca2 = conv2(ts_binarized, pca2_kernel, 'same')';
p.ns_pca2 = conv2(ns_binarized, pca2_kernel, 'same')';

pca3_kernel = interp1(Z.binCenters, score(:, 3), Z.binCenters(1): cfg_master.dt: Z.binCenters(end));
pca3_kernel = pca3_kernel .* gausswin(length(pca3_kernel), 3)';
p.ts_pca3 = conv2(ts_binarized, pca3_kernel, 'same')';
p.ns_pca3 = conv2(ns_binarized, pca3_kernel, 'same')';

% saccade ampl
if cfg_master.skip_amp ~= 1
    ns_idx = find(ns_binarized); ts_idx = find(ts_binarized);
    ns_ampl = ns_binarized; ns_ampl(ns_idx) = nasalAmplitudes;
    ts_ampl = ts_binarized; ts_ampl(ts_idx) = temporalAmplitudes;
    
    ampl_kernel = ones(size(pca1_kernel)); % NOTE: this can end up overlapping with adjacent saccades quite fast...
    p.ns_ampl = conv2(ns_ampl, ampl_kernel, 'same')';
    p.ts_ampl = conv2(ts_ampl, ampl_kernel, 'same')';
end

%% loop over cells
% cellCounterToUse = 0;


fprintf('Cell %d/%d...\n',1,1);
%     cellCounterToUse = cellCounterToUse + 1;
% figure out if we should include this cell

% dependent variable: binned spike train
spk_binned = histc(S.t{1}, TVECe); spk_binned = spk_binned(1:end - 1);
spk_binned = logical(spk_binned > 0); % binarize

p.spk = spk_binned;
p.spk = conv2(p.spk, cfg_master.gausswin, 'same') ./ cfg_master.dt;

% get AHV tc
this_sdf = tsd(sd.TVECc, p.spk');

cfg_tc = [];
cfg_tc.binEdges = cfg_master.tc_binEdges;
cfg_tc.dt = median(diff(sd.TVECc));
tc = TuningCurvesSDF(cfg_tc, this_sdf, this_ahv);
out.tc = tc.tc;

cfg_tc.binEdges = cfg_master.pupil_tc_binEdges;
tc = TuningCurvesSDF(cfg_tc, this_sdf, this_pupilX);
out.pupil_tc = tc.tc;

% saccade predictors based on PETHs
this_xc = xcorr(p.spk, ns_binarized, cfg_master.maxlag, 'unbiased');
out.ns_peth = this_xc;

this_xc = this_xc .* gausswin(length(this_xc), 3); % smooth edges, works better as regressor
ns_conv = conv2(ns_binarized, this_xc', 'same');

this_xc = xcorr(p.spk, ts_binarized, cfg_master.maxlag, 'unbiased');
out.ts_peth = this_xc;

this_xc = this_xc .* gausswin(length(this_xc), 3);
ts_conv = conv2(ts_binarized, this_xc', 'same');

p.ns = ns_conv';
p.ts = ts_conv';
p.sacc = ns_conv' + ts_conv';
out.p = p; 

%this_m = fitglm(p, sd.m.(mn{4}).modelspec, 'Distribution', 'binomial')
for iM = 1:length(mn)
    fprintf('Model %s:\n', mn{iM});
    this_m = fitglm(p, sd.m.(mn{iM}).modelspec, 'Distribution', 'poisson')
    
    out.(mn{iM}).rsq = this_m.Rsquared.Adjusted;
end

out.this_m{iNeuron} = this_m; 
ts_pca1 = this_m.Coefficients{6,1};
ns_pca1 = this_m.Coefficients{7,1};

out.cellname = S.label; % this is the name of the eye tracking neurons (SSN-TT__.t, etc.)

if cfg_master.plotAllPoints == 1
    % get AHV Tuning Curve
    cfg_AHV = [];
    cfg_AHV.subsample_factor = 10;
    [AHV_tsd, tc_out] = AHV_tuning(cfg_AHV, S);
    AHV_dt = median(diff(AHV_tsd.tvec));
    % Firing Rate x AHV
    cfg_Q = [];
    cfg_Q.smooth = 'gauss';
    cfg_Q.gausswin_sd = 0.05;
    cfg_Q.dt = AHV_dt;
    cfg_Q.tvec_edges = AHV_tsd.tvec(1):AHV_dt:AHV_tsd.tvec(end);
    F = MakeQfromS(cfg_Q, S); % convert to FR
    % convert to FR
    F.data = F.data ./ cfg_Q.dt;
    % find FR corresponding to each AHV sample
    F_idx = nearest_idx3(AHV_tsd.tvec, F.tvec);
    AHV_F = F.data(1,F_idx);
    ymax = max(AHV_F);
    out.ahvscatterX = AHV_F;
    out.ahvscatterY = AHV_tsd.data;
end

if cfg_master.debug
    this_fit = this_m.predict(p);
    
    plot(sd.TVECc, p.spk);
    hold on;
    plot(sd.TVECc, this_fit, 'r');
    plot(sd.TVECc(ns_bin), 0, '.c', 'MarkerSize', 20);
    plot(sd.TVECc(ts_bin), 0, '.m', 'MarkerSize', 20);
    set(gca, 'FontSize', 18); ylabel('Firing rate (Hz)'); xlabel('time (s)');
    
%     pause;
end

end