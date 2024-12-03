function [CW, CCW, Z] = AHV_criterion_check(X,Z, tfilelist)
%2024-11-26. JJS.
%   Determines whether each neuron in tfilelist meets the criteria (sinlgy) for r,slope,shuffle(rank), and jointly for all three.
%   Criteria are based on Graham et al., 2023.
% CW.r, CW.slope, CW.rank
% CCW.r, CCW.slope, CCW.rank

CW.r = zeros(1,length(tfilelist));
CW.slope = zeros(1,length(tfilelist));
CW.rank = zeros(1,length(tfilelist));
CW.pass = zeros(1,length(tfilelist));

CCW.r = zeros(1,length(tfilelist));
CCW.slope = zeros(1,length(tfilelist));
CCW.rank = zeros(1,length(tfilelist));
CCW.pass = zeros(1,length(tfilelist));

Z.pass = zeros(1,length(tfilelist));

assert(length(X.rCW) == length(Z.rank_CW) &&  length(X.rCW) == length(tfilelist))

for iNeuron = 1:length(tfilelist)
    if abs(X.bCW(iNeuron)) >= 0.025; CW.slope(iNeuron) = 1; end
    if abs(X.bCCW(iNeuron)) >= 0.025; CCW.slope(iNeuron) = 1; end
    
    if abs(X.rCW(iNeuron)) >= 0.5; CW.r(iNeuron) = 1; end
    if abs(X.rCCW(iNeuron)) >= 0.5;  CCW.r(iNeuron) = 1; end
        
    if Z.CW_pass(iNeuron) == 1; CW.rank(iNeuron) = 1; end
    if Z.CCW_pass(iNeuron) == 1; CCW.rank(iNeuron) = 1; end
    
    if CW.slope(iNeuron) && CW.r(iNeuron) && CW.rank(iNeuron); CW.pass(iNeuron) = 1; end
    if CCW.slope(iNeuron) && CCW.r(iNeuron) && CCW.rank(iNeuron); CCW.pass(iNeuron) = 1; end    
    
    if CW.pass(iNeuron) == 1 || CCW.pass(iNeuron) == 1; Z.pass(iNeuron) = 1; end
end

Z.CW = CW; 
Z.CCW = CCW; 
        
        
        
        