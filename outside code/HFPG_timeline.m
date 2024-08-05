%ylabels = {'Aim 1: behaviour','Aim 1: recording','Aim 1: analysis','Aim 2: behaviour','Aim 2: recording','Aim 2: analysis','Aim 3a (inactivations)','Aim 3b: rec+inact','Aim 3b: analysis'};
ylabels = {'HFPG Aim 1', ...
    'HFPG Aim 2', ...
    'K01 Aim 1', ...
    'K01 Aim 2'};


xlabels = {'Winter 2024','Spring 2025','Fall 2025','Summer 2025','Fall 2025'};

nL = length(ylabels);

rect_ystart = 0.75:1:0.75+nL;
rect_height = 0.5*ones(size(ylabels));
rect_xstart = [0  0  0  1  0  0  1  2];
rect_width =  [3  3  1  2  3  1  1  1];

cols = {'g','g','b','b','b','r','r','r','r','r','r','r'};

figure; hold on;

% draw year lines
for iLine = 1:3
   h = line([iLine iLine],[0.25 nL+0.75]); set(h,'LineStyle','--','Color','k','LineWidth',1); 
end

% make bars
for iL = 1:nL
    
   h = rectangle('Position',[rect_xstart(iL) rect_ystart(iL) rect_width(iL) rect_height(iL)]);
   set(h,'EdgeColor',[0.5 0.5 0.5],'FaceColor',cols{iL});
    
end


set(gca,'YTick',1:nL,'YTickLabel',ylabels,'XTick',0.5:2.5,'XTickLabel',xlabels,'XLim',[0 3],'FontSize',20,'YLim',[0.25 nL+0.75]);
axis ij;

%%
set(gcf,'PaperUnits','inches','PaperSize',[10,4],'PaperPosition',[0 0 10 4])
print('-dpng','-r300','timeline.png')