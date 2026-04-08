function [] = saccadePETHsubplotswZscoreAndAHV(temporaldata, nasaldata, binCenters, TC_all, varargin)

% 4/2021. JJS.
% Make subplots with temporal and nasal saccades overlaid for visual inspection.
xtight = .03;
ytight = .01;
FontSize = 12;
process_varargin(varargin);
numCells = size(temporaldata,1);
numPlots = ceil(sqrt(numCells));
clf
formatSpec = '%.2f';
for iPlot = 1:numCells
    disp(num2str(iPlot))
    t = subtightplot(numPlots,numPlots,iPlot,[xtight ytight]);
%     t = tiledlayout(numPlots,numPlots);
%     ax1 = axes(t);
    hold on
    ax1 = axes;
    a = plot(binCenters, temporaldata(iPlot,:), 'Color', [0    0.4470    0.7410]);
    b = plot(binCenters, nasaldata(iPlot,:), 'Color', [0.8500    0.3250    0.0980]);
%     ax1.XColor = [0.8500    0.3250    0.0980];
%     ax1.YColor = [0.8500    0.3250    0.0980];
    
    ax2 = axes;
    c = plot(TC_all.bins(1,:),smoothdata(TC_all.tc(iPlot,:)), 'Color', [.5 .5 .5]);
%     ax2.XAxisLocation = 'top';
%     ax2.YAxisLocation = 'right';
%     ax2.Color = 'none';
%     ax1.Box = 'off';
%     ax2.Box = 'off';   
    
%     if iPlot == 1
%         legend([a b], 'Temporal', 'Nasal', 'FontSize', 12)
%     end
%     t1 = text(.1,.85,num2str(zVal_T(iPlot),formatSpec), 'Units', 'Normalized', 'Color', a.Color, 'FontSize', FontSize, 'FontWeight', 'bold');
%     t2 = text(.1,.65,num2str(zVal_N(iPlot),formatSpec), 'Units', 'Normalized', 'Color', b.Color, 'FontSize', FontSize);
    title(num2str(iPlot))
end









