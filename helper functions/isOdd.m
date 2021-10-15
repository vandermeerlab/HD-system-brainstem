% isOdd determines if a number is odd or even.
% isOdd takes an array of numbers, and returns a logical array of the same
% size as the input, where true (1) is odd, and false (0) is even.
%
% tf = isOdd(5)
% tf = 
%      1
% tf = isOdd([1 2 3]) 
% tf = 
%      1 0 1
% Created by David Coventry, 8/2/2017

function tf = isOdd(x)

tf = mod(x,2) == 1;

end

