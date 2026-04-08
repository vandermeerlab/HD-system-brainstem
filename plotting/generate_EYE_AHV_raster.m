sd = LoadSessionData([]);

clf
cfg_mr = [];
cfg_mr.SpikeHeight = 1.2;
cfg_mr.lfpWidth = 4;
cfg_mr.lfp(2) = sd.AHV;
diffH = tsd(sd.tsdH.tvec(2:end), diff(sd.tsdH.data)');     % Should this be (1:end-1) or (2:end)?
cfg_mr.lfp(1) = diffH; 
% eye = sd.tsdH;
% eye.data = eye.data';
% cfg_mr.lfp(1) = eye; 
cfg_mr.openNewFig = 0;
h = MultiRasterJeff(cfg_mr, sd.S);

c = axis;
nCells = length(sd.S.t);
axis([c(1) c(2) c(3) c(4)+ nCells])

SSN = HD_GetSSN; disp(SSN);
title(SSN);
set(gca, 'YTick', [])
set(gca, 'YLabel', [])
set(gca, 'FontSize', 25)