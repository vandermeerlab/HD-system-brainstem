% JJS. 9/7/2022.
% For plotting lines from the manually defined stationary periods. 

SSN = HD_GetSSN; disp(SSN);
close all
open(strcat(SSN, '-saccades-edited.fig'))
load(strcat(SSN, '-AHV_stationaryTimes.mat'))
c = axis; 
for iLine = 1:length(STtstart)
    line([STtstart(iLine) STtstart(iLine)], [c(3) c(4)], 'Color', 'c', 'LineWidth', 4, 'LineStyle', '--')  % start of epoch
    line([STtend(iLine) STtend(iLine)], [c(3) c(4)], 'Color', 'm', 'LineWidth', 4, 'LineStyle', '--')      % end of epoch 
    
    line([STtstart(iLine) STtend(iLine)], [c(3) c(4)], 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')      % end of epoch
end

% disp(strcat('num nasal stationary = ', num2str(