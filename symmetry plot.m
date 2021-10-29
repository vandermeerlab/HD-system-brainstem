one = X.bNeg(~Y.F) > 0 & X.bPos(~Y.F) > 0; sone = sum(one);

two = X.bPos(~Y.F) < 0 & X.bNeg(~Y.F) > 0; stwo = sum(two); 

three = X.bNeg(~Y.F) < 0 & X.bPos(~Y.F) < 0; sthree = sum(three);

four = X.bPos(~Y.F) > 0 & X.bNeg(~Y.F) < 0; sfour = sum(four); 
9