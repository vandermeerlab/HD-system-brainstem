function [nasal_MOVING_amp, nasal_REST_amp, temporal_MOVING_amp, temporal_REST_amp, nasal_MOVING_cat, temporal_MOVING_cat, nasal_REST_cat, temporal_REST_cat] = calcSaccadeAmplitudeAllSessions(fd, varargin)

doPlot = 1;
if isempty(fd)
    fd = FindFiles('*keys.m');
end
disp(length(fd))
startSess = 1;
endSess = length(fd);
process_varargin(varargin);

for iSess = startSess: endSess
    pushdir(fileparts(fd{iSess}));
    [nasal_MOVING_amp{iSess}, nasal_REST_amp{iSess}, temporal_MOVING_amp{iSess}, temporal_REST_amp{iSess}] = calcSaccadeAmplitudeSingleSession;
    
end

nasal_MOVING_cat = horzcat(nasal_MOVING_amp{:});
temporal_MOVING_cat = horzcat(temporal_MOVING_amp{:});
nasal_REST_cat = horzcat(nasal_REST_amp{:});
temporal_REST_cat = horzcat(temporal_REST_amp{:});

if doPlot == 1
    bar([median(nasal_MOVING_cat) median(nasal_REST_cat) median(temporal_MOVING_cat) median(temporal_REST_cat)]); hold on
    set(gca, 'XTickLabel', {'nasal MOV'; 'nasal STAT'; 'temporal MOV'; 'temporal STAT'}, 'FontSize', 20)
    ylabel('saccade amplitude', 'FontSize', 20)
    errorbar(1, nanmean(nasal_MOVING_cat), nanstd(nasal_MOVING_cat), 'Color', 'k')
    errorbar(2, nanmean(nasal_REST_cat), nanstd(nasal_REST_cat), 'Color', 'k')
    errorbar(3, nanmean(temporal_MOVING_cat), nanstd(nasal_MOVING_cat), 'Color', 'k')
    errorbar(4, nanmean(temporal_REST_cat), nanstd(nasal_REST_cat), 'Color', 'k')
    title('All sessions')
end