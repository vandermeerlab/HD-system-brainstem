function [] = findrotation()
% JJS. 5/2022.
% Find periods where the platform is stationary by thresholding AHV, and check by plotting. 
threshold = 0.5;

cfg_in = [];
[AHV_tsd] = Get_AHV(cfg_in);

idxS = AHV_tsd.data < threshold & AHV_tsd.data > -threshold;  % tvec indices for Stationary periods. This will include transition points of CW-CCW rotations [need to exclude] 
% St = AHV_tsd.tvec(idxS); 
% dSt = diff(St); 
% dStpoints = dSt > 1; 
% [fr, ~] = find(dStpoints); 
% changepoints = St(fr);
%% exclude periods 



%%
plot(AHV_tsd.tvec, AHV_tsd.data, '.')
hold on
plot(AHV_tsd.tvec(idxS), AHV_tsd.data(idxS), 'r.')