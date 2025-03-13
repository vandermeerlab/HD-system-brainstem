function [Z, numCells] = AHV_criterion_check_ipsi_contra(X, IC, tb, tfilelist)
%2024-11-26. JJS.
%   Determines whether each neuron in tfilelist meets the criteria (sinlgy) for r,slope,shuffle(rank), and jointly for all three.
%   Criteria are based on Graham et al., 2023.
% ipsi.r, ipsi.slope, ipsi.rank
% contra.r, contra.slope, contra.rank
% Input from [CWs,CCWs,slopeToUse,turn_index] = graham_turnbias_index(C, tb, tfilelist)

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
    % This is already calculated in the variables IC.ipsi_pass & IC.contra_pass
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
end
assert(sum(isnan(Z.ipsi_pass3))==0)
assert(sum(isnan(Z.contra_pass3))==0)
assert(sum(isnan(Z.either_pass3))==0)

%% Calculate the fraction of non-AHV, symmetric, asymmetric, asymmetric-unresponsive, and inverted symmetric cells, like Graham et al. Fig. 4C. 

numNonAHV = sum(Z.either_pass3 == 0); perNonAHV = numNonAHV/numCells;

asymmetric = Z.either_pass3 & tb.asymmetric; numAsymmetric = sum(asymmetric); perAsymmetric = numAsymmetric/numCells;

asymmetric_unresponsive = Z.either_pass3 & tb.asymmetric_unresponsive; numAsymmetric_unresponsive = sum(asymmetric_unresponsive); perAsymmetric_unresponsive = numAsymmetric_unresponsive/numCells; 

symmetric = Z.either_pass3 & tb.symmetric_cells; numSymmetric = sum(symmetric); perSymmetric = numSymmetric/numCells; perSymmetric = numSymmetric/numCells;

% there were zero inverted symmetric cells

% bar([perNonAHV perSymmetric perAsymmetric_unresponsive perAsymmetric 0])


labels = {'non-AHV', 'Symmetric', 'Asym-Unresp.', 'Asymmetric', 'Inverted'};

% pie([numNonAHV numSymmetric numAsymmetric_unresponsive numAsymmetric 0])
pie([numNonAHV numSymmetric numAsymmetric_unresponsive numAsymmetric 0], [1 1 1 1 1], labels)



