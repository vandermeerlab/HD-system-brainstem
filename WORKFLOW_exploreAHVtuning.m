%% cd to data folder

%% get AHV
cfg_AHV = [];
cfg_AHV.subsample_factor = 10;
[AHV_tsd, S, tc_out] = AHV_tuning(cfg_AHV);

AHV_dt = median(diff(AHV_tsd.tvec));

%% plot tuning curves

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
AHV_F = F.data(F_idx);

plot(AHV_tsd.data, AHV_F, '.');

%% acf
cfg_acf = [];
cfg_acf.binsize = 0.001;
cfg_acf.max_t = 0.5;
cfg_acf.smooth = 1; % set to 1 to compute ccf on SDF, 0 for raw spikes
cfg_acf.gauss_w = 1; % width of Gaussian convolution window (in s)
cfg_acf.gauss_sd = 0.005; % SD of Gaussian convolution window (in s)
[acf, tvec] = ccf(cfg_acf, S.t{1}, S.t{1});

midpoint = ceil(length(acf)./2);
acf(midpoint) = 0;
plot(tvec, acf);
%%
% cfg_in.dt = .05;
% F = tsd(Q.tvec, Q.data.*1/cfg_in.dt);
% 
% tcmin = min(tc_out.tc);
% tcmax = max(tc_out.tc);
% binsToUse = tcmin:(tcmax-tcmin)/(tc_out.usr.nBins-1):tcmax;
