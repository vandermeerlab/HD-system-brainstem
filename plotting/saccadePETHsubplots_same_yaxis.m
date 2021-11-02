function [] = saccadePETHsubplots_same_yaxis(temporaldata, nasaldata, binCenters, zVal_T, zVal_N, varargin)

% 4/2021. JJS.
% Make subplots with temporal and nasal saccades overlaid for visual inspection.
xtight = .02;
ytight = .01;
FontSize = 14;
process_varargin(varargin);
numCells = size(temporaldata,1);
numPlots = ceil(sqrt(numCells));
clf
formatSpec = '%.2f';
for iPlot = 1:numCells
    disp(num2str(iPlot))
    subtightplot(numPlots,numPlots,iPlot,[xtight ytight])
    set(gca, 'FontSize', 10)
%     yyaxis left
%     set(gca, 'YTick', [])
%     set(gca, 'YTickLabels', [])
%     yticklabels('')
    a = plot(binCenters, temporaldata(iPlot,:), 'r');
    hold on
%     yyaxis right
    b = plot(binCenters, nasaldata(iPlot,:), 'g');
    if iPlot == numCells
        legend([a b], 'Temporal', 'Nasal', 'FontSize', 12)
    end
%     set(gca, 'YTick', [])
    if iPlot < numCells - sqrt(numCells)     
        set(gca, 'XTick', [])
%     else
%         ax = get(gca, 'XTickLabel'); 
%         set(gca, 'XTicklabel', ax, 'FontSize', 10)
    end
%     text(.1,.85,num2str(zVal_T(iPlot),formatSpec), 'Units', 'Normalized', 'Color', a.Color, 'FontSize', FontSize);
%     text(.1,.65,num2str(zVal_N(iPlot),formatSpec), 'Units', 'Normalized', 'Color', b.Color, 'FontSize', FontSize);
    title(num2str(iPlot))
end





