function [] = FindAHVThresholdCrossings(AHV_tsd)

SSN = HD_GetSSN;
% f = FindFiles('*VT1_proc.mat'); 
% load(f{1}, 'pupil');

splitthresh = 100;

idxh = AHV_tsd.data > 1;
idxl = AHV_tsd.data < -1;

idxhtvec = AHV_tsd.tvec(idxh);
idxltvec = AHV_tsd.tvec(idxl); 

IDXH = find(idxh);
diffidxh = diff(IDXH);
tStartH = find(diffidxh>100);
StartPointsH = IDXH(tStartH); 

IDXL = find(idxl);
diffidxl = diff(IDXL);
tStartL = find(diffidxl>splitthresh);
StartPointsL = IDXL(tStartL); 

idxz = AHV_tsd.data < .2 & AHV_tsd.data > -.2;  
IDXz = find(idxz);

clf
plot(AHV_tsd.tvec, AHV_tsd.data); hold on
plot(AHV_tsd.tvec(idxh), AHV_tsd.data(idxh), 'r.')
plot(AHV_tsd.tvec(idxl), AHV_tsd.data(idxl), 'g.')
plot(AHV_tsd.tvec(StartPointsH), zeros(1, length(AHV_tsd.tvec(StartPointsH))), 'h', 'LineWidth', 2, 'Color', 'm')   % magenta stars
plot(AHV_tsd.tvec(StartPointsL), zeros(1, length(AHV_tsd.tvec(StartPointsL))), 'p', 'LineWidth', 2, 'Color', 'b')   % blue stars 
plot(AHV_tsd.tvec(idxz), ones(1, length(AHV_tsd.tvec(idxz))), 'c.')  % 
