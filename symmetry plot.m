one = X.bNeg(~Y.ns) > 0 & X.bPos(~Y.ns) > 0; sone = sum(one);

two = X.bPos(~Y.ns) < 0 & X.bNeg(~Y.ns) > 0; stwo = sum(two); 

three = X.bNeg(~Y.ns) < 0 & X.bPos(~Y.ns) < 0; sthree = sum(three);

four = X.bPos(~Y.ns) > 0 & X.bNeg(~Y.ns) < 0; sfour = sum(four); 

patch([0 1.2 1.2 0],[0 0 1.2 1.2], [.95 .95 .95])
patch([0 -1.2 -1.2 0],[0 0 -1.2 -1.2], [.95 .95 .95])