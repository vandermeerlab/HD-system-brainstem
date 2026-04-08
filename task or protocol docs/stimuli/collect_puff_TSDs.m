function [PUFF_out, CONTROL_out, PUFF_out_sub, CONTROL_out_sub, puff_trials, puff_baseline, control_trials, control_baseline, puff_sub, control_sub, tvec] = collect_puff_TSDs()

bin_for_average = 95; % bin to use for the baseline average of the trial
doPlot = 0;
fd = FindFiles('*.mp4');
cfg_in.plot_time_series = 0;
cfg_in.plot_PETH = 0;
CONTROL_out = [];
PUFF_out = [];
PUFF_out_sub = [];
CONTROL_out_sub = [];
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    [puff_out, control_out, puff_trials, control_trials, cfg_out] = airpuff_peth(cfg_in);
    puff_baseline(iSess) = mean(puff_out.data(1:bin_for_average));
    control_baseline(iSess) = mean(control_out.data(1:75));
    
    PUFF_out = vertcat(PUFF_out, puff_out.data);
    CONTROL_out = vertcat(CONTROL_out, control_out.data);
    
    puff_out_sub = puff_out.data - puff_baseline(iSess);
    PUFF_out_sub = vertcat(PUFF_out_sub, puff_out_sub);
    
    control_out_sub = control_out.data - control_baseline(iSess);
    CONTROL_out_sub = vertcat(CONTROL_out_sub, control_out_sub);
    
    tvec = puff_out.tvec;
end

if doPlot
    clf
    subplot(1,2,1)
    plot(tvec, PUFF_out_sub(1,:)', 'Color', 'r', 'LineWidth', 3)
    hold on
    plot(tvec, CONTROL_out_sub(1,:)', 'Color', 'g', 'LineWidth', 3)
    
    %     plot(tvec, PUFF_out_sub(2,:)', 'Color', 'g')
    %     plot(tvec, PUFF_out_sub(3,:)', 'Color', 'b')
    %     plot(tvec, PUFF_out_sub(4,:)', 'Color', 'm')
    %     legend('Sess 1','Sess 2','Sess 3','Sess 4' )
    set(gca, 'FontSize', 32)
    xlabel('time (sec)')
    ylabel('eye surface area (pixels)')
    
    subplot(1,2,2)   % normalized
    plot(tvec, -PUFF_out_sub(1,:)'./min(PUFF_out_sub(1,:)), 'Color', 'r', 'LineWidth', 3)
    hold on
    plot(tvec, -CONTROL_out_sub(1,:)'./min(CONTROL_out_sub(1,:)), 'Color', 'g', 'LineWidth', 3)
    
    %     plot(tvec, -PUFF_out_sub(2,:)'./min(PUFF_out_sub(2,:)), 'Color', 'g', 'LineWidth', 3)
    %     plot(tvec, -PUFF_out_sub(3,:)'./min(PUFF_out_sub(3,:)), 'Color', 'b', 'LineWidth', 3)
    %     plot(tvec, -PUFF_out_sub(4,:)'./min(PUFF_out_sub(4,:)), 'Color', 'm', 'LineWidth', 3)
%     legend('Sess 1','Sess 2','Sess 3','Sess 4' )
    set(gca, 'FontSize', 32)
    xlabel('time (sec)')
    ylabel('eye surface area (pixels)')
end

puff_sub = [];
control_sub = [];
% trial loop
for iPuff = 1:size(puff_trials.data,1)
    puff_sub(iPuff,1:length(tvec)) = puff_trials.data(iPuff,:) - mean(puff_trials.data(iPuff,1:bin_for_average));
    plot(puff_sub(iPuff,1:length(tvec))); 
    title(num2str(iPuff))
    pause
end

for iControl = 1:size(control_trials.data,1)
    control_sub(iControl,1:length(tvec)) = control_trials.data(iControl,:) - mean(control_trials.data(iControl,1:bin_for_average));
end
