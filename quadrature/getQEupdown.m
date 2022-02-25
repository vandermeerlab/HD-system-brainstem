function updownTSD = getQEupdown(chA,chB)
% JJS. 2022-02-25. 
% getQUupdown.m  Takes wheel encoder input--two CSC channels--and gives a tsd output of ones and zeros, indicating 'high' or 'low' voltage values, respectively. 
% Output is a single 2-D tsd. This function assumes that channel A (chA) is the first channel of the quadrature encoder and channel B (chB) is the second. 
% Order is important for interpreting direction. 

%   Inputs
%       chA: first csc from the wheel encoder  
%       chB: second csc from the wheel encoder 
%   Outputs
%       updownTSD:  tsd of ones and zeros which indicates up and down states, respectively. row 1 = chA. row 2 = chB. tvec is CSC timestamps. 

if nargin ~=2
     disp('Loading default CSC channels') 
     cfg = []; 
     filenameA = FindFile('*CSC34.Ncs');     % for sessions with wheel encoder data, CSC34 is the first output channel of the encoder.
     [~, name, ext] = fileparts(filenameA);
     cfg.fc = {strcat(name,ext)}; 
     chA = LoadCSC(cfg); 
     
     cfg = [];
     filenameB = FindFile('*CSC35.Ncs');     % for sessions with wheel encoder data, CSC35 is the first output channel of the encoder.
     [~, name, ext] = fileparts(filenameB); 
     cfg.fc = {strcat(name,ext)}; 
     chB = LoadCSC(cfg); 
end

assert(length(chA.tvec) == length(chB.tvec))   % size of the tsd's should be identical
% Error testing: NaNs
if sum(isnan(chA.data))>0; warning('NaNs detected in chA. Shouldnt be any NaNs'); end
if sum(isnan(chB.data))>0; warning('NaNs detected in chB. Shouldnt be any NaNs'); end

meanA = mean(chA.data);   % find the mean. Up states should be well above the mean. Down states should be well below the mean. 
meanB = mean(chB.data); 

% Error testing: Values equal to the mean
numA = sum(chA.data==meanA); if numA > 0; warning('one or more values in chA exactly equal to the mean voltage'); end
numB = sum(chB.data==meanB); if numB > 0; warning('one or more values in chA exactly equal to the mean voltage'); end

vectorA = chA.data > meanA;    % logical with value of 1's for putative 'up' states and 0's for putative 'down' states
vectorB = chB.data > meanB;

updownData = vertcat(vectorA, vectorB); 
updownTSD = tsd(chA.tvec, updownData);












