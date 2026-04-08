function [ipsi_all_three_pass, contra_all_three_pass, AHV_responsive, AHV_table, ipsiversive, contraversive] = check_three_criteria(X, IC)
%2024-12-02. JJS. Check all three of Jalina's criteria to identify a neuron as AHV-sensitive. These criteria are slope, r, and shuffle test (of r values). 
% Relies on output of AHV_pearson_correlation2.m 
%   Inputs
%               X - structure with key fields like slopeIsGood, rIsGood, and bothIsGood. 
numCells = length(X.slopeToUse);

for iNeuron = 1:numCells
    if X.IPSI_slope_and_r_IsGood(iNeuron) & IC.ipsi_pass(iNeuron)
        ipsi_all_three_pass(iNeuron) = 1;
    else
        ipsi_all_three_pass(iNeuron) = 0;
    end
    
    if X.CONTRA_slope_and_r_IsGood(iNeuron) & IC.contra_pass(iNeuron) %#ok<*AND2>
        contra_all_three_pass(iNeuron) = 1;
    else
        contra_all_three_pass(iNeuron) = 0;
    end
    
    if ipsi_all_three_pass(iNeuron) | contra_all_three_pass(iNeuron)
        AHV_responsive(iNeuron) = 1;
    else
        AHV_responsive(iNeuron) = 0;
    end   
end
ipsi_slope_greater = X.slopeToUse == 1;
contra_slope_greater = X.slopeToUse == 0;

ipsiversive = AHV_responsive & ipsi_slope_greater; 
contraversive = AHV_responsive & contra_slope_greater; 

assert(sum(ipsiversive) + sum(contraversive) == sum(AHV_responsive));

AHV_table.numAHVresponsive = sum(AHV_responsive);
AHV_table.per_ipsi_pass = sum(ipsi_all_three_pass)/numCells;
AHV_table.per_contra_pass = sum(contra_all_three_pass)/numCells;
AHV_table.per_any_pass = sum(AHV_responsive)/numCells;

AHV_table.num_ipsiversive = sum(ipsiversive);
AHV_table.num_contraversive = sum(contraversive);

AHV_table.per_ipsiversive = AHV_table.num_ipsiversive/AHV_table.numAHVresponsive;
AHV_table.per_contraversive = AHV_table.num_contraversive/AHV_table.numAHVresponsive;



