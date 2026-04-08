function [w] = calc_slopes_r_sparse(tc_out)
%2024-11-23. JJS.
%   Calculate the slope, r value, rsq, and p value for each side of the AHV tuning curve.
x = tc_out.binCenters; x = x';  % tc_out.binCenters(34) is the innermost negative bin [-6 to 0]. tc_out.binCenters(35) is the innermost positive bin [0 to +6].
y = tc_out.tc;    y = y';
idx = isnan(y);
x(:,2) = ones(length(x),1);
%% Corr for Positive AHV values (CCW)
CCWindex = tc_out.binCenters > 0 & tc_out.binCenters <= 90;
[bCCW,~,~,~,statsCCW] = regress(y(CCWindex),x(CCWindex,:));  % b is the slope (spk/s/deg/s). stats is Rsq, F, p-value, & error variance. Rsq is the coefficient of determination.
w.slopeCCW = bCCW(1);
w.CCWyInt = bCCW(2);
w.RsqCCW = statsCCW(1);
w.pCCW = statsCCW(3);
temp2 = corrcoef(x(CCWindex' & ~idx,1),y(CCWindex' & ~idx));
w.rCCW = temp2(1,2);   % this is the r value (correlation coefficient)

%% Corr for Negative AHV values (CW)
CWindex = tc_out.binCenters < 0 & tc_out.binCenters >=-90;
[bCW,~,~,~,statsCW] = regress(y(CWindex),x(CWindex,:));
w.slopeCW = bCW(1);
w.CWyInt = bCW(2);
w.RsqCW = statsCW(1);
w.pCW = statsCW(3);
temp3 = corrcoef(x(CWindex' & ~idx,1),y(CWindex' & ~idx));
w.rCW = temp3(1,2);

end