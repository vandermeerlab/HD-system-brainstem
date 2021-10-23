

% Inputs:   X  -  structre with fields...
%                       rsqAll:     [1×n double]     Rsq statistic for each neuron, fit with a single line
%                       pAll:       [1×n double]     p-values for each neuron, fit with a single line
%                       coeffsAll:  [n×2 double]     1st row = slope (Beta value); 2nd row = y-intercept. Single line fit 
%                       rsqPos:     [1×n double]     Rsq for positive part of tuning curve (above zero deg/sec)
%                       pPos:       [1×n double]     p-values for positive part of tuning curve 
%                       coeffsPos:  [n×2 double]     coeffs for positive part of tuning curve
%                       rsqNeg:     [1×n double]     Rsq for negative part of tuning curve
%                       pNeg:       [1×n double]     p-values for negative part of tuning curve
%                       coeffsNeg:  [n×1 double]     coeffs for negative part of tuning curve 


indexP = X.pPos<thrsehold & X.pNeg<thrsehold;

indexSamePos = X.coeffsPos(:,1)>=0 &  X.coeffsNeg(:,1)>=0; 

indexSameNeg = X.coeffsPos(:,1)<=0 &  X.coeffsNeg(:,1)<=0; 


symmetric = indexSamePos | indexSameNeg; 

asymmetric1 = indexSamePos & indexSameNeg; 

symmetricSig = symmetric & indexP'; 