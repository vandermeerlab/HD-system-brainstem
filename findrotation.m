function [] = findrotation()
% JJS. 12/2021.
% Find periods where the platform is stationary by thresholding AHV, and check by plotting. 
minPoints = 100; 

cfg_in = [];
[AHV_tsd] = Get_AHV(cfg_in);

idxS = AHV_tsd.data < .2 & AHV_tsd.data > -.2;  % Stationary periods 
St = AHV_tsd.tvec(idxS); 
dSt = diff(St); 
dStpoints = dSt > 10; 
[fr, fc] = find(dStpoints); 

idxR = ~idxS; 