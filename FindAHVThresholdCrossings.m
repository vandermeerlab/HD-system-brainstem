splitthresh = 100;

idxh = AHVtsd.data > 1;
idxl = AHVtsd.data < -1;

idxhtvec = AHVtsd.tvec(idxh);
idxltvec = AHVtsd.tvec(idxl); 

IDXH = find(idxh);
diffidxh = diff(IDXH);
tStartH = find(diffidxh>100);
StartPointsH = IDXH(tStartH); 
gitgit s

IDXL = find(idxl);
diffidxl = diff(IDXL);
tStartL = find(diffidxl>splitthresh);
StartPointsL = IDXL(tStartL); 


clf
plot(AHVtsd.tvec, AHVtsd.data); hold on
plot(AHVtsd.tvec(idxh), AHVtsd.data(idxh), 'r.')
plot(AHVtsd.tvec(idxl), AHVtsd.data(idxl), 'g.')
plot(AHVtsd.tvec(StartPointsH), zeros(1, length(AHVtsd.tvec(StartPointsH))), 'h', 'LineWidth', 2, 'Color', 'm')   % magenta stars
plot(AHVtsd.tvec(StartPointsL), zeros(1, length(AHVtsd.tvec(StartPointsL))), 'p', 'LineWidth', 2, 'Color', 'b')   % blue stars 