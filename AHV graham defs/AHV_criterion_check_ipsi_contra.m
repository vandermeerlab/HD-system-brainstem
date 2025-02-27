function [Z, numCells] = AHV_criterion_check_ipsi_contra(X, IC, tfilelist)
%2024-11-26. JJS.
%   Determines whether each neuron in tfilelist meets the criteria (sinlgy) for r,slope,shuffle(rank), and jointly for all three.
%   Criteria are based on Graham et al., 2023.
% ipsi.r, ipsi.slope, ipsi.rank
% contra.r, contra.slope, contra.rank

numCells = length(tfilelist);

Z.ipsi.r_pass = zeros(1,length(tfilelist));
Z.ipsi.slope_pass = zeros(1,length(tfilelist));
Z.ipsi.rank_pass = zeros(1,length(tfilelist));
Z.ipsi.pass = zeros(1,length(tfilelist));

Z.contra.r_pass = zeros(1,length(tfilelist));
Z.contra.slope_pass = zeros(1,length(tfilelist));
Z.contra.rank_pass = zeros(1,length(tfilelist));
Z.contra.pass = zeros(1,length(tfilelist));

Z.ipsi_pass3 = NaN(1, length(tfilelist));
Z.contra_pass3 = NaN(1, length(tfilelist));
Z.either_pass3 = NaN(1, length(tfilelist));

for iNeuron = 1:length(tfilelist)
    % Slope
    if abs(X.ipsi.slope(iNeuron)) >= 0.025; Z.ipsi.slope_pass(iNeuron) = 1; end
    if abs(X.contra.slope(iNeuron)) >= 0.025; Z.contra.slope_pass(iNeuron) = 1; end
    % r value
    if abs(X.ipsi.r(iNeuron)) >= 0.5; Z.ipsi.r_pass(iNeuron) = 1; end
    if abs(X.contra.r(iNeuron)) >= 0.5;  Z.contra.r_pass(iNeuron) = 1; end
    % Rank
    
    %% All 3 - ipsi
    if Z.ipsi.slope_pass(iNeuron) && Z.ipsi.r_pass(iNeuron) && IC.ipsi_pass(iNeuron)
        Z.ipsi_pass3(iNeuron) = 1;
    else
        Z.ipsi_pass3(iNeuron) = 0;
    end
    
    %% All 3 - contra
    if Z.contra.slope_pass(iNeuron) && Z.contra.r_pass(iNeuron) && IC.contra_pass(iNeuron)
        Z.contra_pass3(iNeuron) = 1;
    else 
        Z.contra_pass3(iNeuron) = 0;
    end
    
    %% All 3 - either
    if Z.ipsi_pass3(iNeuron) == 1 || Z.contra_pass3(iNeuron) == 1
        Z.either_pass3(iNeuron) = 1; 
    else 
        Z.either_pass3(iNeuron) = 0; 
    end
    
    %% Identify Symmetric Cells
    
    
    
    
    
    
end
assert(sum(isnan(Z.ipsi_pass3))==0)
assert(sum(isnan(Z.contra_pass3))==0)
assert(sum(isnan(Z.either_pass3))==0)

