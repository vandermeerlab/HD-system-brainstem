function [data_out, data_outM, data_outD, data_outMD] = getSlowPhaseData2(cfg_in)
% function [data_out, nSpikesRemoved] = getSlowPhaseData(cfg_in, sd, iCell)
%
%   Inputs
%           sd          - init file with sd.S
%           iCell       - which cell number to plot. leave empty [] if not calculating neural data
%           cfg_in      - config variable. Choose the "de-saccading" window.
%   Outputs
%           data_out    - the TSDs for this session, including horiz_eye_pos, horiz_eye_vel, horiz_eye_vel_smooth, AHV, and HD (a.k.a. "orientation")
%           data_outD   - the same TSDs, but with the pre and post saccade interval cut out (with antirestrict.m)
%           data_outMD  - same as above, but additionally with the stationary periods removed (with antirestrict.m)
%           
%                           data_out.HD                         - tsd of HD (platform orientation), in common timebase
%                           data_outM.HD                        - tsd of HD (platform orientation), without stationary periods
%                           data_outD.HD                        - tsd of HD (platform orientation), without saccade epochs
%                           data_outMD.HD                       - tsd of HD (platform orientation), without stationary periods or saccade epochs

%                           data_out.AHV                        - tsd of AHV, in common timebase
%                           data_outM.AHV                       - tsd of AHV, without stationary epochs
%                           data_outD.AHV                       - tsd of AHV, without saccade epochs
%                           data_outMD.AHV                      - tsd of AHV, without stationary periods or saccade epochs

%                           data_out.horiz_eye_pos              - tsd of EP, in common timebase
%                           data_outM.horiz_eye_pos             - tsd of EP, without stationary periods
%                           data_outD.horiz_eye_pos             - tsd of EP, without saccade epochs
%                           data_outMD.horiz_eye_pos            - tsd of EP, without stationary periods or saccade epochs

%                           data_out.horiz_eye_vel              - tsd of EV, in common timebase 
%                           data_outM.horiz_eye_vel             - tsd of EV, without stationary periods 
%                           data_outD.horiz_eye_vel             - tsd of EV, without saccade epochs 
%                           data_outMD.horiz_eye_vel            - tsd of EV, without stationary periods or saccade epochs

%                           data_out.horiz_eye_vel_smooth       - tsd of filtered & smoothed EV, in common timebase 
%                           data_outM.horiz_eye_vel_smooth      - tsd of filtered & smoothed EV, without stationary periods
%                           data_outD.horiz_eye_vel_smooth      - tsd of filtered & smoothed EV, without saccade epochs 
%                           data_outMD.horiz_eye_vel_smooth     - tsd of filtered & smoothed EV, without stationary periods or saccade epochs

%                           data_out.wheel                      - tsd of wheel data, in common timebase
%                           data_outM.wheel                     - tsd of wheel data, without stationary periods
%                           data_outD.wheel                     - tsd of wheel data, without saccade epochs
%                           data_outMD.wheel                    - tsd of wheel data, without stationary periods or saccade epochs


% Uses saccade .mat file to restrict data to slow phase only
% Useful for later tuning curve and GLM analysis
%
% 2024-04-18. MvdM, based on JJS getSlowPhaseTC3
% 2024-05-29. JJS. Added conditions to skip the resampling of neural data.
% 2024-09-04. JJS. Updated to use the variables that are saved in units of degrees and degrees/second (tsdHdeg, diffHdeg).
% 2024-09-08. JJS. This version (ver2) no longer takes sd or iCell as inputs. sd can take a long time to load, and there is no need to do anything with neurons at this step.
start = tic; 
data_out = []; data_outD = []; data_outMD = [];

cfg_def = [];
cfg_def.FontSize = 20;
cfg_def.saccade_pre = 0.1;  % how many seconds to cut out before the saccade.
cfg_def.saccade_post = 0.1; % how many seconds to cut out after the saccade.
cfg_def.doPlot = 1;
cfg_def.markerSize = 3;
cfg_def.plotWhichPhase = 'slow'; % 'slow', 'both'
cfg_def.medfilt_window_size = 11; % number of samples to use for slow phase velocity filter
cfg_def.smooth = 1;
cfg_def.doWarn = 1;
cfg_def.downsamplefactor = 10;
cfg_def.verbose = 1;

cfg_Q = [];
cfg_Q.dt = 0.02; % binsize in s
cfg_Q.smooth = 'gauss'; % [], 'gauss'
cfg_Q.gausswin_size = 1; % gaussian window size in seconds; only used if cfg.smooth = 'gauss'
cfg_Q.gausswin_sd = 0.1; % SD for gaussian convolution in seconds
cfg_def.cfg_Q = cfg_Q;

cfg_master = ProcessConfig2(cfg_def, cfg_in);

%% Setup
SSN = HD_GetSSN; disp(SSN)
load(strcat(SSN, '-saccades-edited.mat'), 'AHV_tsd', 'combinedSaccades', 'tsdHdeg', 'diffHdeg') % load the edited (curated) saccade timestamps and the pupil trace time series (position and velocity)
load(strcat(SSN, '-AHV_StationaryTimes.mat'), 'STtstart', 'STtend')           % load the manually generated stationary period start/stop times

timestouse = ~isnan(combinedSaccades);                                      % create a logical of non-NaN saccade times. ***At some point, go through and remove the NaNs.
combinedSaccadesToUse = combinedSaccades(timestouse);                       % select only non-NaN saccade values to use

pre = combinedSaccadesToUse - cfg_master.saccade_pre;                       % get pre-saccade timestamps
post = combinedSaccadesToUse + cfg_master.saccade_post;                     % get post-saccade timestamps

if cfg_master.doWarn == 1; warning('off','all'); end                        % remove warnings, especially having to do with resampleTSD
tsdHdeg.data = tsdHdeg.data';                                               % .data needs to be 1 x n format to be a "well-formed tsd" for resampleTSD.m
diffHdeg.tvec = diffHdeg.tvec';                                             % .tvec needs to be n x 1 format to be a "well-formed tsd". *** Need to resave these vars so they have the correct (and consistent) format

%% RESAMPLE data to be on common timebase
cfg.method = 'rebinning';
new_tvec_edges = AHV_tsd.tvec(1): cfg_master.cfg_Q.dt: AHV_tsd.tvec(end); % take AHV_tsd from saccades-edited.mat file. This will start at t0 = 0 and tEnd = .tvec(end)
new_tvec = new_tvec_edges(1:end-1) + (cfg_master.cfg_Q.dt / 2);
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%       WHEEL
[wheel] = get_wheel_speed;
tic; 
if cfg_master.verbose == 1; disp('rebinning WHEEL'); end 
wheel_resamp = resampleTSD(cfg, wheel, new_tvec);
toc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%       HEAD DIRECTION
tic; 
if cfg_master.verbose == 1; disp('rebinning HD'); end
[~, orientation, ~, ~, ~, ~] = GetOrientationValues([]);
[orientationtouse] = downsampleOrientationValues(orientation, cfg_master.downsamplefactor); % if downsample factor = 10, it will go from 2kHz to 200Hz
head_direction_resamp = resampleTSD(cfg, orientationtouse, new_tvec);
toc; 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%        AHV
tic; 
if cfg_master.verbose == 1; disp('rebinning AHV'); end
cfg = []; cfg.method = 'rebinning';
ahv_resamp = resampleTSD(cfg, AHV_tsd, new_tvec);
toc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%       EYE POSITION
tic; 
if cfg_master.verbose == 1; disp('rebinning EP'); end
horiz_eye_pos_resamp = resampleTSD(cfg, tsdHdeg, new_tvec);
toc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%       EYE VELOCITY
tic; 
if cfg_master.verbose == 1; disp('rebinning EV'); end
velocity_filt = medfilt1(diffHdeg.data, cfg_master.medfilt_window_size);                % apply median filter
velocity_smooth = smoothdata(velocity_filt);                                            % apply smoothing
horiz_eye_vel_smooth = tsd(diffHdeg.tvec, velocity_smooth);                             % make tsd
horiz_eye_vel_resamp = resampleTSD(cfg, diffHdeg, new_tvec);
horiz_eye_vel_smooth_resamp = resampleTSD(cfg, horiz_eye_vel_smooth, new_tvec);
toc

% Save resampled data
data_out.wheel = wheel_resamp;                                                          % WHEEL
data_out.HD = head_direction_resamp;                                                    % HD
data_out.AHV = ahv_resamp;                                                              % AHV
data_out.horiz_eye_pos = horiz_eye_pos_resamp;                                          % EP 
data_out.horiz_eye_vel = horiz_eye_vel_resamp;                                          % EV
data_out.horiz_eye_vel_smooth = horiz_eye_vel_smooth_resamp;                            % EV, smoothed & filtered

%% Remove stationary periods only (M = "moving")
wheelM                  = antirestrict(wheel_resamp, STtstart, STtend);                 % WHEEL 
head_directionM         = antirestrict(head_direction_resamp, STtstart, STtend);        % HD
ahvM                    = antirestrict(ahv_resamp, STtstart, STtend);                   % AHV
horiz_eye_posM          = antirestrict(horiz_eye_pos_resamp, STtstart, STtend);         % EP
horiz_eye_velM          = antirestrict(horiz_eye_vel_resamp, STtstart, STtend);         % EV
horiz_eye_vel_smoothM   = antirestrict(horiz_eye_vel_smooth_resamp, STtstart, STtend);  % EV, smoothed & filtered

% Save the resampled, MOVING-only data 
data_outM.wheel                 = wheelM;                                               % WHEEL
data_outM.HD                    = head_directionM;                                      % HD
data_outM.AHV                   = ahvM;                                                 % AHV
data_outM.horiz_eye_pos         = horiz_eye_posM;                                       % EP
data_outM.horiz_eye_vel         = horiz_eye_velM;                                       % EV
data_outM.horiz_eye_vel_smooth  = horiz_eye_vel_smoothM;                                % EV, smooth & filtered 

%% Remove quick phase periods only (D = "de-saccaded")
wheelR                  = antirestrict(wheel_resamp, pre, post);                        % WHEEL 
head_directionR         = antirestrict(head_direction_resamp, pre, post);               % HD
ahvR                    = antirestrict(ahv_resamp, pre, post);                          % AHV
horiz_eye_posR          = antirestrict(horiz_eye_pos_resamp, pre, post);                % EP
horiz_eye_velR          = antirestrict(horiz_eye_vel_resamp, pre, post);                % EV
horiz_eye_vel_smoothR   = antirestrict(horiz_eye_vel_smooth_resamp, pre, post);         % EV, smoothed & filtered

% Save the resampled, DESACCADED data
data_outD.wheel                 = wheelR;                                               % WHEEL
data_outD.HD                    = head_directionR;                                      % HD
data_outD.AHV                   = ahvR;                                                 % AHV
data_outD.horiz_eye_pos         = horiz_eye_posR;                                       % EP
data_outD.horiz_eye_vel         = horiz_eye_velR;                                       % EV
data_outD.horiz_eye_vel_smooth  = horiz_eye_vel_smoothR;                                % EV, smooth & filtered 

%% Remove quick phase periods from the MOVING-only data (M = "moving" + D = "de-saccaded")
wheelMD                  = antirestrict(wheelM, pre, post);                        % WHEEL 
head_directionMD         = antirestrict(head_directionM, pre, post);               % HD
ahvMD                    = antirestrict(ahvM, pre, post);                          % AHV
horiz_eye_posMD          = antirestrict(horiz_eye_posM, pre, post);                % EP
horiz_eye_velMD          = antirestrict(horiz_eye_velM, pre, post);                % EV
horiz_eye_vel_smoothMD   = antirestrict(horiz_eye_vel_smoothM, pre, post);         % EV, smoothed & filtered

% Save the resampled, DESACCADED data
data_outMD.wheel                 = wheelMD;                                               % WHEEL
data_outMD.HD                    = head_directionMD;                                      % HD
data_outMD.AHV                   = ahvMD;                                                 % AHV
data_outMD.horiz_eye_pos         = horiz_eye_posMD;                                       % EP
data_outMD.horiz_eye_vel         = horiz_eye_velMD;                                       % EV
data_outMD.horiz_eye_vel_smooth  = horiz_eye_vel_smoothMD;                                % EV, smooth & filtered 

%% Plot 
if cfg_master.doPlot
    figure(1); clf;
    displayFactor = 10;
    plot(tsdHdeg.tvec, tsdHdeg.data .* displayFactor); hold on  % pupil position
    plot(diffHdeg.tvec, diffHdeg.data, 'r');  % pupil velocity
    plot(horiz_eye_vel_smooth.tvec, horiz_eye_vel_smooth.data, 'g', 'LineWidth', 2) % smoothed pupil velocity
    set(gca, 'FontSize', cfg_master.FontSize)
    plot(tsdHdeg.tvec, tsdHdeg.data .* displayFactor, '.b') % pupil position as points
    plot(diffHdeg.tvec, diffHdeg.data, 'r.') % pupil velocity as points
    plot(horiz_eye_vel_smooth.tvec, horiz_eye_vel_smooth.data, 'g', 'LineWidth', 5) % smoothed pupil velocity. Replotted here so the green line is in front.
    legend('pupil position', 'velocity', 'median filtered & smoothed velocity')
    title(SSN)
    xlabel('time (s)')
    ylabel('deg/s')
end

if cfg_master.doWarn == 1; warning('on','all'); end  % Return warnings to ON state. 
toc(start)
