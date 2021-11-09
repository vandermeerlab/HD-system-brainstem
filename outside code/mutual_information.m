function I = mutual_information(X,Y,px,py,pxy)
% Simple function to pull out the mutual information between x and y in
% bits.
% Assumes that sum_i p(x_i)==1
%              sum_j p(y_j)==1
%              sum_j sum_i p(x_i,y_j)==1
%
%

% I(X;Y) = sum_{y E Y} sum_{x E X} p(x,y) log_2(p(x,y)/(p(x)p(y))

uniqueX = unique(X);
uniqueY = unique(Y);
[uniqueXY,id] = unique([X Y],'rows');
px = px(id);
py = py(id);
pxy = pxy(id);
X = X(id);
Y = Y(id);

% Assume px sums to 1 over all unique x.
prob = zeros(length(uniqueX),1);
for ix = 1 : length(uniqueX)
    idx = uniqueX(ix)==X;
    prob(ix) = mean(px(idx));
end
prob = prob/sum(prob);
for ix = 1 : length(uniqueX)
    idx = uniqueX(ix)==X;
    px(idx) = prob(ix);
end

% Assume py sums to 1 over all unique y.
prob = zeros(length(uniqueY),1);
for iy = 1 : length(uniqueY)
    idy = uniqueY(iy)==Y;
    prob(iy) = mean(py(idy));
end
prob = prob/sum(prob);
for iy = 1 : length(uniqueY)
    idy = uniqueY(iy)==Y;
    py(idy) = prob(iy);
end

% Assume pxy sums to 1 over all unique XY pairs.
prob = zeros(size(uniqueXY,1),1);
for ixy = 1 : size(uniqueXY,1)
    idxy = uniqueXY(ixy,1)==X&uniqueXY(ixy,2)==Y;
    prob(idxy) = mean(pxy(idxy));
end
prob = prob/sum(prob);
for ixy = 1 : size(uniqueXY,1)
    idxy = uniqueXY(ixy,1)==X&uniqueXY(ixy,2)==Y;
    pxy(idxy) = prob(ixy);
end


P = pxy.*log(pxy./(px.*py))./log(2);
I = sum(P);