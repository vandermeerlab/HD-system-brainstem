SSN = HD_GetSSN;
f = FindFiles('*VT1_proc.mat'); 
load(f{1}, 'pupil');


splitthresh = 100;

idxh = AHVtsd.data > 1;
idxl = AHVtsd.data < -1;

idxhtvec = AHVtsd.tvec(idxh);
idxltvec = AHVtsd.tvec(idxl); 

IDXH = find(idxh);
diffidxh = diff(IDXH);
tStartH = find(diffidxh>100);
StartPointsH = IDXH(tStartH); 

IDXL = find(idxl);
diffidxl = diff(IDXL);
tStartL = find(diffidxl>splitthresh);
StartPointsL = IDXL(tStartL); 

idxz = AHVtsd.data < .2 & AHVtsd.data > -.2;  
IDXz = find(idxz);

clf
plot(AHVtsd.tvec, AHVtsd.data); hold on
plot(AHVtsd.tvec(idxh), AHVtsd.data(idxh), 'r.')
plot(AHVtsd.tvec(idxl), AHVtsd.data(idxl), 'g.')
plot(AHVtsd.tvec(StartPointsH), zeros(1, length(AHVtsd.tvec(StartPointsH))), 'h', 'LineWidth', 2, 'Color', 'm')   % magenta stars
plot(AHVtsd.tvec(StartPointsL), zeros(1, length(AHVtsd.tvec(StartPointsL))), 'p', 'LineWidth', 2, 'Color', 'b')   % blue stars 
plot(AHVtsd.tvec(idxz), ones(1, length(AHVtsd.tvec(idxz))), 'c.')  % 
