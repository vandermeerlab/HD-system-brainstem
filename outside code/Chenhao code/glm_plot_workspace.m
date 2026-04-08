%%
% eye sensitivity 2D plot
% ax = figure('Position', get(0, 'Screensize'));
FontSize = 20;
clf
for i = 1:N_glmNPH
    disp(num2str(i))
    if Rsq_bestNPH(i) <= 0 || isnan(Rsq_bestNPH(i))
        Rsq_bestNPH(i) = 1e-6;
    end
    scatter(EP_sensitivity_bestNPH(i,1),EV_sensitivity_bestNPH(i,1),Rsq_bestNPH(i)*1000,...
        'MarkerFaceColor','black',...
        'MarkerFaceAlpha',0,...
        'MarkerEdgeColor','black',...
        'MarkerEdgeAlpha',1);
    hold on;
    errorbar(EP_sensitivity_bestNPH(i,1),EV_sensitivity_bestNPH(i,1),...
        EP_sensitivity_bestNPH(i,2),'horizontal','black');
    errorbar(EP_sensitivity_bestNPH(i,1),EV_sensitivity_bestNPH(i,1),...
        EV_sensitivity_bestNPH(i,2),'vertical','black');
end
xlabel('EP sensitivity ((sp/s)/deg)');
ylabel('EV sensitivity ((sp/s)/(deg/s))');

scatter(5,-1.4,1*1000,...
    'MarkerFaceColor','black',...
    'MarkerFaceAlpha',0,...
    'MarkerEdgeColor','black',...
    'MarkerEdgeAlpha',1);
text(5.5,-1.4,'R2=1', 'FontSize', FontSize);

scatter(5,-1.8,0.5*1000,...
    'MarkerFaceColor','black',...
    'MarkerFaceAlpha',0,...
    'MarkerEdgeColor','black',...
    'MarkerEdgeAlpha',1);
text(5.5,-1.8,'R2=0.5', 'FontSize', FontSize);

scatter(5,-2.2,0.25*1000,...
    'MarkerFaceColor','black',...
    'MarkerFaceAlpha',0,...
    'MarkerEdgeColor','black',...
    'MarkerEdgeAlpha',1);
text(5.5,-2.2,'R2=0.25', 'FontSize', FontSize);

rectangle('Position',[4.3, -2.4, 3.2, 1.25]);

% c = axis; 
% line([-10 10], [-10 10], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
set(gca, 'FontSize', 36)

savefig(['glm_fig_35EVthreshold_scatter.fig']);

%%
% 1D plot
% ax = figure('Position', get(0, 'Screensize'));
LineWidth = 3;
clf
subplot(141);
bar(0,mean(bias_bestNPH(:,1)),10);
hold on;
e = errorbar(0,mean(bias_bestNPH(:,1)),std(bias_bestNPH(:,1)),'black', 'LineWidth', 3);
h = scatter(randn(length(bias_bestNPH(:,1)),1),bias_bestNPH(:,1), 'filled', 'LineWidth', 5)
set(gca, 'LineWidth', LineWidth)
xticks([]);
xlim([-10,10]);
ylim([0,250]);
xlabel('Bias');
ylabel('spk/s');
set(gca, 'FontSize', FontSize)

subplot(142);
bar(0,mean(abs(EP_sensitivity_bestNPH(:,1))),10);
hold on;
errorbar(0,mean(abs(EP_sensitivity_bestNPH(:,1))),std(abs(EP_sensitivity_bestNPH(:,1))),'black', 'LineWidth', 3);
scatter(randn(length(EP_sensitivity_bestNPH(:,1)),1),abs(EP_sensitivity_bestNPH(:,1)), 'filled', 'LineWidth', 5);
set(gca, 'LineWidth', LineWidth)
xticks([]);
xlim([-10,10]);
ylim([0,10]);
xlabel('EP sensitivity');
ylabel('sp/s/deg');
set(gca, 'FontSize', FontSize)
c = axis;
axis([c(1) c(2) c(3) 7.6]);

subplot(143);
bar(0,mean(abs(EV_sensitivity_bestNPH(:,1))),10);
hold on;
errorbar(0,mean(abs(EV_sensitivity_bestNPH(:,1))),std(abs(EV_sensitivity_bestNPH(:,1))),'black', 'LineWidth', 3);
scatter(randn(length(EV_sensitivity_bestNPH(:,1)),1),abs(EV_sensitivity_bestNPH(:,1)), 'filled', 'LineWidth', 5);
set(gca, 'LineWidth', LineWidth)
xticks([]);
xlim([-10,10]);
ylim([0,4]);
xlabel('EV sensitivity');
ylabel('sp/s/deg/s');
set(gca, 'FontSize', FontSize)
c = axis;
axis([c(1) c(2) c(3) 7.6]);

subplot(144);
bar(0,mean(abs(AHV_sensitivity_bestNPH(:,1))),10);
hold on;
errorbar(0,mean(abs(AHV_sensitivity_bestNPH(:,1))),std(abs(AHV_sensitivity_bestNPH(:,1))),'black', 'LineWidth', 3);
scatter(randn(length(AHV_sensitivity_bestNPH(:,1)),1),abs(AHV_sensitivity_bestNPH(:,1)), 'filled', 'LineWidth', 5);
set(gca, 'LineWidth', LineWidth)
xticks([]);
xlim([-10,10]);
ylim([0,4]);
xlabel('AHV sensitivity');
ylabel('sp/s/deg/s');
set(gca, 'FontSize', FontSize)
c = axis;
axis([c(1) c(2) c(3) 7.6]);









savefig(['glm_fig_35EVthreshold_EyeSensitivity.fig']);

%% 
% FR_residue with AHV
ax = figure();
bar(0,mean(abs(AHV_sensitivity_bestNPH(:,1))),10);
hold on;
errorbar(0,mean(abs(AHV_sensitivity_bestNPH(:,1))),std(abs(AHV_sensitivity_bestNPH(:,1))),'black');
scatter(randn(length(AHV_sensitivity_bestNPH(:,1)),1),abs(AHV_sensitivity_bestNPH(:,1)),'k.');
xticks([]);
xlim([-10,10]);
ylim([0,5]);
xlabel('AHV sensitivity');
ylabel('AHV sensitivity ((sp/s)/(deg/s))');

savefig(['glm_fig_35EVthreshold_AHVSensitivity.fig']);

%%
% R2
% ax = figure();
clf
bar(10,mean(abs(AHV_Rsq_bestNPH(:,1))),10);
hold on;
bar(-10,mean(abs(Rsq_bestNPH(:,1))),10);
errorbar(10,mean(abs(AHV_Rsq_bestNPH(:,1))),std(abs(AHV_Rsq_bestNPH(:,1))),'black');
errorbar(-10,mean(abs(Rsq_bestNPH(:,1))),std(abs(Rsq_bestNPH(:,1))),'black');
scatter(10,abs(AHV_Rsq_bestNPH(:,1)),'k.');
scatter(-10,abs(Rsq_bestNPH(:,1)),'k.');
for i = 1:length(Rsq_bestNPH(:,1))
    plot([-10,10],[Rsq_bestNPH(i,1),AHV_Rsq_bestNPH(i,1)],'k');
end
% xticks([-10,10]);
xlim([-20,20]);
ylim([0,1]);
% xticklabels({'p<.01^-20','FR residue ~ 1 + AHV residue, conditioning on EP and EV'});
set(gca, 'XTicklabels', [])

ylabel('R2');

savefig(ax, ['glm_fig_35EVthreshold_R2.fig']);