function [out,dt] = FindNextEventAfter(ref,findset)
% function [out,dt] = FindNextEventAfter(ref,findset)

l = length(ref);
out = nan(l,1);
dt = nan(l,1);

for iE = 1:l
   
    temp = findset-ref(iE); temp(temp < 0) = NaN;
    [val,ind] = min(temp);
    out(iE) = findset(ind); 
    dt(iE) = out(iE)-ref(iE);
    
end