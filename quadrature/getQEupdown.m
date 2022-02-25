function updownTSD = getQEupdown(cfg_in)
% Takes wheel encoder input--two CSC channels--and gives a tsd output of ones and zeros, indicating 'high' or 'low' voltage values, respectively. 
% Output is a single 2-D tsd. This function assumes that channel A (chA) is the first channel of the quadrature encoder and channel B (chB) is the second. 
% Order is important for interpreting direction. 
%
% JJS. 2022-02-25. 
% MvdM edit

%   Outputs:
%       updownTSD:  tsd of ones and zeros which indicates up and down states, respectively. row 1 = chA. row 2 = chB. tvec is CSC timestamps. 
%       use as input to XXXX
%   cfg_in defaults
%       cfg_in.chA = '*CSC34.Ncs'; % filename of quadrature encoder channel 1
%       cfg_in.chB = '*CSC35.Ncs'; % filename of quadrature encoder channel 2

cfg_def = [];
cfg_def.verbose = 1; 
cfg_in.chA = '*CSC34.Ncs'; % filename of quadrature encoder channel 1
cfg_in.chB = '*CSC35.Ncs'; % filename of quadrature encoder channel 2

cfg = ProcessConfig(cfg_def, cfg_in);

%
please = []; 
please.fc{1} = FindFile(cfg.chA);
chA = LoadCSC(please); 

please.fc{1} = FindFile(cfg.chB);
chB = LoadCSC(please); 
     
assert(length(chA.tvec) == length(chB.tvec))   % size of the tsd's should be identical

meanA = mean(chA.data); 
meanB = mean(chB.data); 

vectorA = chA.data > meanA;    % logical with value of 1's for putative 'up' states and 0's for putative 'down' states
vectorB = chB.data > meanB;

updownData = vertcat(vectorA, vectorB); 
updownTSD = tsd(chA.tvec, updownData);

