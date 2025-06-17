function [CWs,CCWs,slopeToUse,turn_index, tb, AHV_preferred_side] = graham_turnbias_index(C, tfilelist)
%2025-02-27. JJS
%   Calculates the "normalized turn bias score" (what I'm calling "asymmetry index"), as outlined in Graham et al., 2023.
%   Input comes from ------------   [X, C, neuronList, numSess, cfg_out] = AHV_pearson_correlation_CONTRA_IPSI(cfg_in, tfilelist)   ---------------------
% 
%        AHV_preferred_side(iNeuron) = -1; % negative AHV is the preferred AHV side. CW cell
%        AHV_preferred_side(iNeuron) = 1; % positive AHV is the preferred AHV side. CCW cell
        
PlotSlopes= 0;        
CWs = C.slopeCW;
CCWs = C.slopeCCW;
Smax = NaN(1,length(tfilelist));

%% Which side has the greater slope?      1 = CW. 0 = CCW
for iNeuron = 1:length(tfilelist)
    if abs(CWs(iNeuron)) > abs(CCWs(iNeuron))
        slopeToUse(iNeuron) = 1;
    elseif abs(CWs(iNeuron)) < abs(CCWs(iNeuron))
        slopeToUse(iNeuron) = 0;
    else
        error('slopes are identical')
    end
end

Smax(slopeToUse==1) = CWs(slopeToUse==1);
Smax(slopeToUse==0) = CCWs(slopeToUse==0);

assert(sum(isnan(Smax))==0);

% turnbias = abs((abs(CWs + CCWs)/2.*Smax));  % this is from Jalina's paper. It does not make sense to me.
% for iNeuron = 1: length(tfilelist)
%     turnbias(iNeuron) = abs(abs(CWs(iNeuron) + CCWs(iNeuron)) / 2 * Smax(iNeuron));  
% end

for iNeuron = 1:length(tfilelist)
    turn_index(iNeuron) = (CWs(iNeuron) + CCWs(iNeuron))/(abs(CWs(iNeuron)) + abs(CCWs(iNeuron)));
end

tb.symmetric_cells = turn_index <0.3 & turn_index > -0.3;
disp(strcat('num symmetric cells = ', num2str(sum(tb.symmetric_cells))))

tb.CW_asymmetric = turn_index <= -0.7;
disp(strcat('num CW asymmetric cells = ', num2str(sum(tb.CW_asymmetric))))

tb.CCW_asymmetric = turn_index >= 0.7;
disp(strcat('num CCW asymmetric cells = ', num2str(sum(tb.CCW_asymmetric))))

tb.CW_asymmetric_unresponsive = turn_index >-0.7 & turn_index <= -0.3;
disp(strcat('num CW asymmetric-unresponsive cells = ', num2str(sum(tb.CW_asymmetric_unresponsive))))

tb.CCW_asymmetric_unresponsive = turn_index >= 0.3 & turn_index < 0.7;
disp(strcat('num CCW asymmetric-unresponsive cells = ', num2str(sum(tb.CCW_asymmetric_unresponsive))))

tb.asymmetric = tb.CCW_asymmetric | tb.CW_asymmetric;
tb.asymmetric_unresponsive = tb.CW_asymmetric_unresponsive | tb.CCW_asymmetric_unresponsive;

assert(sum(tb.symmetric_cells) + sum(tb.CW_asymmetric) + sum(tb.CCW_asymmetric) + sum(tb.CW_asymmetric_unresponsive) + ...
sum(tb.CCW_asymmetric_unresponsive) == length(tfilelist))

tb.num_symmetric = sum(tb.symmetric_cells); 
tb.num_asymmetric = sum(tb.CW_asymmetric) + sum(tb.CCW_asymmetric);
tb.num_asymmetric_unresponsive = sum(tb.CW_asymmetric_unresponsive) + sum(tb.CCW_asymmetric_unresponsive);

%% Label the AHV side with higher firing
for iNeuron = 1:length(turn_index)
    if turn_index(iNeuron) < 0
        AHV_preferred_side(iNeuron) = -1; % negative AHV is preferred side. CW cell
        
        Preferred_slope(iNeuron) = CWs(iNeuron);
        nonPreferred_slope(iNeuron) = CCWs(iNeuron);
    elseif turn_index(iNeuron) >= 0 
        AHV_preferred_side(iNeuron) = 1; % positive AHV is preferred side. CCW cell
        
        Preferred_slope(iNeuron) = CCWs(iNeuron);
        nonPreferred_slope(iNeuron) = CWs(iNeuron);
    else
        error('issue with turn index')
    end
end

if PlotSlopes
    plot(Preferred_slope, nonPreferred_slope, '.')
    
    
       









