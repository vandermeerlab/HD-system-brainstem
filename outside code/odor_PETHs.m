%% multi-unit activity (sum of all spikes over time)
% assumes data has been loaded, as in plot_HC_depth_profile_2.m
cfg_MUA = []; 
%cfg_MUA.tvec = lfp_tsd.tvec'; % timebase to compute MUA on
cfg_MUA.tvec = (firstSpike(S):0.001:lastSpike(S))'; % timebase to compute MUA on
MUA = getMUA(cfg_MUA, S);

figure; % Create a figure
subplot(121) %Create a subplot of channel 121
plot(MUA) %Plot MUA
title('raw MUA (spks/s)') %Add a title 

subplot(122) %Create a subplot of channel 122
MUAz = zscore_tsd(MUA); %Create a variable that produces the zscore
plot(MUAz) %Plot the zscore
title('z-scored MUA (a.u.)'); %Add a title 


%% event-triggered MUA
for iC = 1:3
    cue_evt(iC) = SelectTS([], evt, iC); % select events for each individual cue
end

cfg_peth = []; % parameters for PETH
cfg_peth.window = [-1 1];
cfg_peth.dt = 0.01;
cfg_peth.mode = 'interp';

figure;

cols = 'rgb';
for iC = 1:3

    % remove illegal events
    tstart = cue_evt(iC).t{1} + cfg_peth.window(1); tend = cue_evt(iC).t{1} + cfg_peth.window(2);
    remove_idx = tstart < MUAz.tvec(1) | tend > MUAz.tvec(end);
    cue_evt(iC).t{1}(remove_idx) = [];

    out = TSDpeth(cfg_peth, MUAz, cue_evt(iC));
    h(iC) = plot(out, 'Color', cols(iC), 'LineWidth', 2);
    hold on;
end

set(gca, 'FontSize', 18, 'TickDir', 'out'); box off;
legend(h, {'Cue 1', 'Cue 2', 'Cue 3'}); legend boxoff;

xlabel('Time from cue onset (s)')
ylabel('z-scored MUA')

vline(0,':');

%% some single-unit PETHs
nCells = length(S.t);
for iC = 1:nCells

    this_S = SelectTS([], S, iC); % select only the current cell
    MUA = getMUA(cfg_MUA, this_S); % "MUA" for one cell is just that cell's firing rate
    MUAz = zscore_tsd(MUA);

    figure(1 + floor((iC-1)/9));
    subplot(3,3,mod(iC-1, 9) + 1);

    for iCue = 1:3
        this_out = TSDpeth(cfg_peth, MUAz, cue_evt(iCue)); hold on;
        plot(this_out, 'Color', cols(iCue), 'LineWidth', 2);
        title(sprintf('cell %d', iC));
    end
    vline(0,':');
end