function [data_out, data_outR] = getSlowPhaseData2(cfg_in)
% function [data_out, nSpikesRemoved] = getSlowPhaseData(cfg_in, sd, iCell)
%
%   Inputs
%           sd          - init file with sd.S 
%           iCell       - which cell number to plot. leave empty [] if not calculating neural data
%           cfg_in      - config variable. Choose the "de-saccading" window. 
%   Outputs
%           data_out    - the TSDs for this session, including horiz_eye_pos, horiz_eye_vel, horiz_eye_vel_smooth, AHV, and HD (a.k.a. "orientation")       
%           data_outR   - the same TSDs, but with the pre and post saccade interval cut out (with antirestrict.m) 

% Uses saccade .mat file to restrict data to slow phase only
% Useful for later tuning curve and GLM analysis
%
% 2024-04-18. MvdM, based on JJS getSlowPhaseTC3
% 2024-05-29. JJS. Added conditions to skip the resampling of neural data.
% 2024-09-04. JJS. Updated to use the variables that are saved in units of degrees and degrees/second (tsdHdeg, diffHdeg). 
% 2024-09-08. JJS. This version (ver2) no longer takes sd or iCell as inputs. sd can take a long time to load, and there is no need to do anything with neurons at this step. 

cfg_def = [];
cfg_def.FontSize = 20;
cfg_def.saccade_pre = 0.1;  % how many seconds to cut out before the saccade.
cfg_def.saccade_post = 0.1; % how many seconds to cut out after the saccade.
cfg_def.doPlot = 1;
cfg_def.markerSize = 3;
cfg_def.plotWhichPhase = 'slow'; % 'slow', 'both'
cfg_def.medfilt_window_size = 11; % number of samples to use for slow phase velocity filter
cfg_def.smooth = 1;

cfg_Q = [];
cfg_Q.dt = 0.02; % binsize in s
cfg_Q.smooth = 'gauss'; % [], 'gauss'
cfg_Q.gausswin_size = 1; % gaussian window size in seconds; only used if cfg.smooth = 'gauss'
cfg_Q.gausswin_sd = 0.1; % SD for gaussian convolution in seconds
cfg_def.cfg_Q = cfg_Q;

cfg_master = ProcessConfig2(cfg_def, cfg_in);

SSN = HD_GetSSN; disp(SSN)
if exist(strcat(SSN, '-saccades-edited.mat'),'file') == 2
    load(strcat(SSN, '-saccades-edited.mat')) % this will load the edited (curated) saccade timestamps and the pupil trace time series (position and velocity)
else
    error('saccade mat file not found')
end

%% resample data to be on common timebase
new_tvec_edges = AHV_tsd.tvec(1): cfg_master.cfg_Q.dt: AHV_tsd.tvec(end); % take AHV_tsd from saccades-edited.mat file. 
new_tvec = new_tvec_edges(1:end-1) + (cfg_master.cfg_Q.dt / 2);

% ahv
cfg = []; cfg.method = 'rebinning';
ahv_resamp = resampleTSD(cfg, AHV_tsd, new_tvec);

tsdHdeg.data = tsdHdeg.data'; % .data needs to be 1 x n format to be a "well-formed tsd" for resampleTSD.m
diffHdeg.tvec = diffHdeg.tvec'; % .tvec needs to be n x 1 format to be a "well-formed tsd". *** Need to resave these vars so they have the correct (and consistent) format
 
% eye position
horiz_eye_pos_resamp = resampleTSD(cfg, tsdHdeg, new_tvec); 

% eye velocity
%horiz_eye_vel = dxdt(sd.tsdH.tvec, sd.tsdH.data, 'window', 0.2, 'postSmoothing', 0.1); % this would be nice, why doesn't it work?
velocity_filt = medfilt1(diffHdeg.data, cfg_master.medfilt_window_size);
horiz_eye_vel_smooth = tsd(diffHdeg.tvec, velocity_filt);

if cfg_master.doPlot
    clf;
    displayFactor = 10;
    plot(tsdHdeg.tvec, tsdHdeg.data .* displayFactor); hold on  % pupil position
    plot(diffHdeg.tvec, diffHdeg.data, 'r');  % pupil velocity
    plot(horiz_eye_vel_smooth.tvec, horiz_eye_vel_smooth.data, 'g', 'LineWidth', 2) % smoothed pupil velocity
    set(gca, 'FontSize', cfg_master.FontSize)
    plot(tsdHdeg.tvec, tsdHdeg.data .* displayFactor, '.b') % pupil position as points
    plot(diffHdeg.tvec, diffHdeg.data, 'r.') % pupil velocity as points
    plot(horiz_eye_vel_smooth.tvec, horiz_eye_vel_smooth.data, 'g', 'LineWidth', 2) % smoothed pupil velocity. Replotted here so the green line is in front.
    legend('pupil position', 'velocity', 'smoothed velocity')
end

cfg = []; cfg.method = 'rebinning';
horiz_eye_vel_resamp = resampleTSD(cfg, diffHdeg, new_tvec);
horiz_eye_vel_smooth_resamp = resampleTSD(cfg, horiz_eye_vel_smooth, new_tvec);

% head direction
tic
[~, orientation, ~, ~, ~, ~] = GetOrientationValues([]);
toc
cfg = []; cfg.method = 'rebinning';
head_direction_resamp = resampleTSD(cfg, orientation, new_tvec);

%% Limit the data to everything but the quick phase periods
timestouse = ~isnan(combinedSaccades); % create a logical with not NaN saccade times
combinedSaccadesToUse = combinedSaccades(timestouse); % select only non-NaN saccade values to use

pre = combinedSaccadesToUse - cfg_master.saccade_pre;
post = combinedSaccadesToUse + cfg_master.saccade_post;

%% restrict everything else and set up data_out
ahvR = antirestrict(ahv_resamp, pre, post);
horiz_eye_posR = antirestrict(horiz_eye_pos_resamp, pre, post);
horiz_eye_velR = antirestrict(horiz_eye_vel_resamp, pre, post);
horiz_eye_vel_smoothR = antirestrict(horiz_eye_vel_smooth_resamp, pre, post);
head_directionR = antirestrict(head_direction_resamp, pre, post);

% restricted and resampled data
data_outR.horiz_eye_pos = horiz_eye_posR; data_out.horiz_eye_pos = horiz_eye_pos_resamp;
data_outR.horiz_eye_vel = horiz_eye_velR; data_out.horiz_eye_vel = horiz_eye_vel_resamp;
data_outR.horiz_eye_vel_smooth = horiz_eye_vel_smoothR; data_out.horiz_eye_vel_smooth = horiz_eye_vel_smooth_resamp;
data_outR.AHV = ahvR; data_out.AHV = ahv_resamp;
data_outR.HD = head_directionR; data_out.HD = head_direction_resamp;
