%%
cd('C:\Users\mvdm\Dropbox\projects\Jeff\ahv_peth_pca');
load SaccadePETHs

%% prep data
x1 = FRxBinNsmooth';
x2 = FRxBinTsmooth';
x = cat(2, x1, x2);

nBins = size(x, 1);
nCells = size(x1, 2);

x0 = mean(x);  
%xc = x - repmat(x0, [n 1]);
%xc = xc';
%% svd
%[u, s, v] = svd(xc); 
%evalue = diag(s(1:3,1:3)); 

%figure;
%plot(binCenters, u(:,1));
%hold on;
%plot(binCenters, u(:,2), 'r');
%plot(binCenters, u(:,3), 'g');

%%
[coeff, score, lat, tsq, var_expl, mu] = pca(x);

%% plot some diagnostics
figure;

subplot(221);
plot(binCenters, score(:,1), 'b');
hold on;
plot(binCenters, score(:,2), 'r');
plot(binCenters, score(:,3), 'g');
xlabel('time');
title('first 3 PCs');
set(gca, 'FontSize', 18);

subplot(222);
plot(cumsum(var_expl)); set(gca, 'XLim', [1 10]);
hold on;
plot(cumsum(var_expl), '.b', 'MarkerSize', 20 ); set(gca, 'XLim', [1 10]);
plot([1 10], [95 95], 'r:');
xlabel('PC#');
ylabel('variance explained');
set(gca, 'FontSize', 18);

%%
rec = score(:, 1:50)*coeff(:, 1:50)' + repmat(mu, [nBins 1]);

%% reconstruct specific cell TCs
subplot(223); 

this_cell = 7; 
plot(binCenters, x1(:, this_cell), 'k', 'LineWidth', 2); hold on;
%c0 = [0.8 0.8 0.8];
c = 'brgcm';

this_rec = repmat(mu(this_cell), [nBins 1]);
plot(binCenters, this_rec, 'Color', c0);
for iPC = 1:3
    this_rec = this_rec + score(:, iPC)*coeff(this_cell, iPC)';
    %plot(binCenters, this_rec, 'Color', c0 - iPC.*[0.15 0.15 0.15]);
    plot(binCenters, this_rec, 'Color', c(iPC));
end

xlabel('time');
title(sprintf('c%dn, cf %.2f %.2f %.2f', this_cell, coeff(this_cell,1), coeff(this_cell,2), coeff(this_cell,3)));
set(gca, 'FontSize', 18);

subplot(224);
plot(binCenters, x2(:, this_cell), 'k', 'LineWidth', 2); hold on;
%c0 = [0.8 0.8 0.8];
c = 'brgcm';

this_cell = this_cell + nCells;
this_rec = repmat(mu(this_cell), [nBins 1]);
plot(binCenters, this_rec, 'Color', c0);
for iPC = 1:3
    this_rec = this_rec + score(:, iPC)*coeff(this_cell, iPC)';
    %plot(binCenters, this_rec, 'Color', c0 - iPC.*[0.15 0.15 0.15]);
    plot(binCenters, this_rec, 'Color', c(iPC));
end

xlabel('time');
title(sprintf('c%dt, cf %.2f %.2f %.2f', this_cell - nCells, coeff(this_cell,1), coeff(this_cell,2), coeff(this_cell,3)));
set(gca, 'FontSize', 18);

%%
figure;

subplot(221);

n_idx = 1:nCells;
plot3(coeff(n_idx, 1), coeff(n_idx, 2), coeff(n_idx, 3), '.r', 'MarkerSize', 10);
s1h = gca;
grid on; 
set(gca, 'FontSize', 18);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); title ('nasal');
ax = axis;

subplot(222);

h = plot3(0, 0, 0); delete(h);
s2h = gca;
axis(ax)
for iC = 1:nCells
   h = text(coeff(iC, 1), coeff(iC, 2), coeff(iC, 3), num2str(iC)); set(h, 'FontSize', 10);
   hold on;
end
grid on; 
set(gca, 'FontSize', 18);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); title ('nasal');

Link = linkprop([s1h, s2h],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);


subplot(223);

n_idx = nCells+1:size(coeff, 1);
plot3(coeff(n_idx, 1), coeff(n_idx, 2), coeff(n_idx, 3), '.g', 'MarkerSize', 10);
s3h = gca;
grid on; 
set(gca, 'FontSize', 18);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); title ('temporal');
ax = axis;

subplot(224);

h = plot3(0, 0, 0); delete(h);
s4h = gca;
axis(ax)
for iC = nCells+1:size(coeff, 1)
   h = text(coeff(iC, 1), coeff(iC, 2), coeff(iC, 3), num2str(iC-nCells)); set(h, 'FontSize', 10);
   hold on;
end
grid on; 
set(gca, 'FontSize', 18);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); title ('temporal');

Link = linkprop([s3h, s4h],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);

%% highlight some selected TCs
toPlot = [7 23 31 30 58 47]; % nasal high PC1
toPlot = [56 13 11 32 53 41 43 65]; % temporal high PC1
toPlot = [11 26];
toPlot = [55 48 57 4 18 20 1 36];

figure;
for iP = 1:length(toPlot)
    
    subplot(3, 3, iP);
    this_c = toPlot(iP);
    plot(binCenters, x1(:, this_c), 'r', 'LineWidth', 2); hold on;
    plot(binCenters, x2(:, this_c), 'g', 'LineWidth', 2); hold on;
    title(sprintf('c%d', this_c));
    set(gca, 'FontSize', 16); grid on;
    
end

%%
figure;
subplot(221)

plot(coeff(1:nCells, 1), coeff(nCells + 1:end, 1), '.k', 'MarkerSize', 10);
grid on;
set(gca, 'FontSize', 18);
xlabel('nasal PC1 coeff'); ylabel('temporal PC1 coeff');
s1h = gca; ax = axis;


subplot(222)


h = plot(0, 0); delete(h);
s2h = gca;
axis(ax)
for iC = 1:nCells
   h = text(coeff(iC, 1), coeff(iC+nCells, 1), num2str(iC)); set(h, 'FontSize', 10);
   hold on;
end
grid on;
set(gca, 'FontSize', 18);
xlabel('nasal PC1 coeff'); ylabel('temporal PC1 coeff');