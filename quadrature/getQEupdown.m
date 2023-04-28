function up_down_tsd = getQEupdown(cfg_in)
% Takes wheel encoder input--two CSC channels--and gives a tsd output of ones and zeros, indicating 'high' or 'low' voltage values, respectively. 
% Output is a single 2-D tsd. This function assumes that channel A (chA) is the first channel of the quadrature encoder and channel B (chB) is the second. 
% Order is important for interpreting direction. 
%
% JJS. 2022-02-25. 
% MvdM edit

%   Outputs:
%       up_down_tsd:  tsd of ones and zeros which indicates up and down states, respectively. row 1 = chA. row 2 = chB. tvec is CSC timestamps. 
%       use as input to XXXX
%   cfg_in defaults
%       cfg_in.chA = '*CSC34.Ncs'; % filename of quadrature encoder channel 1
%       cfg_in.chB = '*CSC35.Ncs'; % filename of quadrature encoder channel 2

cfg_def = [];
cfg_def.verbose = 1; 
cfg_def.chA = '*CSC34.Ncs'; % filename of quadrature encoder channel 1   % changed from cfg_in.chA = ... to cfg_def.chA = ...
cfg_def.chB = '*CSC35.Ncs'; % filename of quadrature encoder channel 2

cfg = ProcessConfig(cfg_def, cfg_in);

%
please = []; 
please.fc{1} = FindFile(cfg.chA);
chA = LoadCSC(please); 

please.fc{1} = FindFile(cfg.chB);
chB = LoadCSC(please); 
     
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
up_down_tsd = tsd(chA.tvec, updownData);












