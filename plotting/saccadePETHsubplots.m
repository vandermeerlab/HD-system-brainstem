function [] = saccadePETHsubplots(temporaldata, nasaldata, binCenters, cellID, X, varargin)

% 4/2021. JJS.
% Make subplots with temporal and nasal saccades overlaid for visual inspection. 
xtight = .02;
ytight = .01;
process_varargin(varargin);
numCells = size(temporaldata,1);
numPlots = ceil(sqrt(numCells));
clf
for iPlot = 1:numCells
    disp(num2str(iPlot))
    subtightplot(numPlots,numPlots,iPlot,[xtight ytight])
    yyaxis left
    set(gca, 'YTick', [])
    yticklabels('')
    a = plot(binCenters, temporaldata(iPlot,:));
    c1 = axis;
%     line([0 0], [0 c1(4)], 'Color', 'k', 'LineWidth', 1)
    yyaxis right 
    b = plot(binCenters, nasaldata(iPlot,:));
    if iPlot == numCells
        legend([a b], 'Temporal', 'Nasal', 'FontSize', 12)
    end
    c2 = axis;
%     line([0 0], [0 c2(2)], 'Color', 'k', 'LineWidth', 1)
    set(gca, 'XTick', [])
    set(gca, 'YTick', [])
    
    formatSpec = '%.2f';
    ylim=get(gca,'ylim');
    xlim=get(gca,'xlim');
    text(.9*xlim(1),.9*ylim(2),num2str(X.rsqAll(cellID(iPlot)),formatSpec));
    title(num2str(iPlot))
end





