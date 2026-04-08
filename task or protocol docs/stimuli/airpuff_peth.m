function [puff_out, control_out, puff_trials, control_trials, cfg_out] = airpuff_peth(cfg_in)

cfg_def.plot_time_series = 0;
cfg_def.plot_PETH = 0;

cfg_out = ProcessConfig(cfg_def, cfg_in);


SSN = HD_GetSSN; disp(SSN)
load(strcat(SSN,'-VT1_proc.mat'))
events = LoadEvents([]);
%% Find session start time
index = strfind(events.label, 'Starting Recording');
if index{1} == 4                                 % Start Recording should be in the first or second .label position.
    starttime = events.t{4}(1);  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
elseif index{5} == 1
    starttime = events.t{5}(1); % for session with laser events, the start recording time moves to position 2.
else
    error('could not find start time for this session')
end
%% Load events and subtract start time
Events = [];
for iE = 1:length(events.t)
    Events.t{iE} = events.t{iE} - starttime;  % start record should be event.t{5}
end
Events.type = events.type;
Events.label = events.label;
%% Load video timestamps and subtract start time
[~, b, c] = fileparts(FindFile('*VT1.smi'));
fn = strcat(b,c);
tic
tvec_raw = read_smi(fn);
toc
pupiltime = tvec_raw - starttime;
%% Make easy names for key events
puff_tone_starts = Events.t{2};
puff_tone_ends = Events.t{8};
puff_deliveries = Events.t{1};

control_tone_starts = Events.t{4};
control_tone_ends = Events.t{9};
control_deliveries = Events.t{3};

if cfg_out.plot_time_series
    %% Plot blink time series
    LineWidth = 5;
    clf; plot(pupiltime, blink_0, 'LineWidth', LineWidth); hold on; c = axis;
    set(gca, 'FontSize', 32)
    xlabel('time (sec)')
    ylabel('eye surface area (pixels)')
    line([puff_deliveries(1) puff_deliveries(1)], [c(3) c(4)], 'Color', 'r', 'LineWidth', 5, 'LineStyle', '-') % puff times
    line([control_deliveries(1) control_deliveries(1)], [c(3) c(4)], 'Color', 'g', 'LineWidth', 5, 'LineStyle', '-') % control times
    
    line([puff_deliveries puff_deliveries], [c(3) c(4)], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-') % puff times
    line([control_deliveries control_deliveries], [c(3) c(4)], 'Color', 'g', 'LineWidth', 2, 'LineStyle', '-') % control times
    legend('Blink trace', 'Airpuffs', 'Control puffs')
    title(SSN)
    disp('press any key to continue')
    pause
end

%% Calculate PETHs
cfg_peth = [];
cfg_peth.mode = 'interp';  % b/c M539-2024-09-25-1 session crashes otherwise 
cfg_peth.dt = 0.02;
blink_tsd = tsd(pupiltime, blink_0);
[puff_out, puff_trials]  = TSDpeth(cfg_peth, blink_tsd, puff_deliveries);
[control_out, control_trials]  = TSDpeth(cfg_peth, blink_tsd, control_deliveries);

% mot_tsd = tsd(pupiltime, motion_1);
% [mot_out, mot_trials]  = TSDpeth(cfg_peth, mot_tsd, puff_deliveries);


%% Plot peths
if cfg_out.plot_PETH
    figure
    plot(puff_out.tvec, puff_out.data, 'LineWidth', 5, 'Color', 'r'); hold on
    plot(control_out.tvec, control_out.data, 'LineWidth', 5, 'Color', 'g')
    xlabel('time (sec)')
    ylabel('eye surface area (pixels)')
    set(gca, 'FontSize', 32)
    legend('Airpuff response', 'Control response')
    title(SSN)
end
%% See if there is learning across the session.




