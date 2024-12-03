function [ipsi, contra, IC] = AHV_criterion_check_ipsi_contra(X, IC, tfilelist)
%2024-11-26. JJS.
%   Determines whether each neuron in tfilelist meets the criteria (sinlgy) for r,slope,shuffle(rank), and jointly for all three.
%   Criteria are based on Graham et al., 2023.
% ipsi.r, ipsi.slope, ipsi.rank
% contra.r, contra.slope, contra.rank

ipsi.r_pass = zeros(1,length(tfilelist));
ipsi.slope_pass = zeros(1,length(tfilelist));
ipsi.rank_pass = zeros(1,length(tfilelist));

ipsi.pass = zeros(1,length(tfilelist));

contra.r_pass = zeros(1,length(tfilelist));
contra.slope_pass = zeros(1,length(tfilelist));
contra.rank_pass = zeros(1,length(tfilelist));

contra.pass = zeros(1,length(tfilelist));

IC.pass = zeros(1,length(tfilelist));

% assert(length(X.ipsi.r) == length(IC.r_ipsi_true) &&  length(X.ipsi.r) == length(tfilelist))

for iNeuron = 1:length(tfilelist)
    % Slope 
    if abs(X.ipsi.slope(iNeuron)) >= 0.025; ipsi.slope_pass(iNeuron) = 1; end
    if abs(X.contra.slope(iNeuron)) >= 0.025; contra.slope_pass(iNeuron) = 1; end
    % r value
    if abs(X.ipsi.r(iNeuron)) >= 0.5; ipsi.r_pass(iNeuron) = 1; end
    if abs(X.ipsi.r(iNeuron)) >= 0.5;  contra.r_pass(iNeuron) = 1; end
    % Rank     
    if IC.rank_ipsi(iNeuron) == 1; ipsi.rank_pass(iNeuron) = 1; end
    if IC.rank_contra(iNeuron) == 1; contra.rank_pass(iNeuron) = 1; end
    % All 3 - ipsi
    if ipsi.slope_pass(iNeuron) && ipsi.r_pass(iNeuron) && ipsi.rank_pass(iNeuron); ipsi.pass(iNeuron) = 1; end
    % All 3 - contra 
    if contra.slope_pass(iNeuron) && contra.r_pass(iNeuron) && contra.rank_pass(iNeuron); contra.pass(iNeuron) = 1; end    
    % All 3 - either
    if ipsi.pass(iNeuron) == 1 || contra.pass(iNeuron) == 1; IC.pass(iNeuron) = 1; end
end

IC.ipsi = ipsi; 
IC.contra = contra; 
        
        
        
        