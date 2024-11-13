function [out] = SESSION_ahv_glmfit_tfilelist_AHV_only(cfg_in, tfilelist)
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

[path, neuron_to_use, ext] = fileparts(tfilelist);
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

%%
p = table; % contains regressors

nMaxVars = 5; % only used for initializing t-stat matrix
% % baseline model MUST be defined first or things will break!
sd.m.baseline.modelspec = 'spk ~ 1 + time';
% %sd.m.ahv.modelspec = 'spk ~ 1 + time + ns';

sd.m.ahv_true.modelspec = 'spk ~ 1 + time + ahv';  % pos. and neg. AHV values
 
sd.m.ahv_absolute.modelspec = 'spk ~ 1 + time + ahv_abs';  % absolute value of AHV [for looking at SYMMETRIC tuning]

sd.m.ahv_both.modelspec = 'spk ~ 1 + time + ahv + ahv_abs';  % both 'true' and 'absolute' tuning together 

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

p.ahv_abs = interp1(AHV_tsd.tvec, abs(AHV_tsd.data), sd.TVECc)';

fprintf('Cell %d/%d...\n',1,1);

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

out.p = p;

%this_m = fitglm(p, sd.m.(mn{4}).modelspec, 'Distribution', 'binomial')
for iM = 1:length(mn)
    fprintf('Model %s:\n', mn{iM});
    [this_m] = fitglm(p, sd.m.(mn{iM}).modelspec, 'Distribution', 'poisson') %#ok<NOPRT>
    
    out.(mn{iM}).rsq = this_m.Rsquared.Adjusted;
end

out.this_m = this_m;
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
