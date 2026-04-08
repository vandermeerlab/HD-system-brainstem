function [CWs,CCWs,slopeToUse] = graham_asymmetric_index(X, C, tfilelist)
%2025-02-27. JJS 
%   Calculates the "normalized turn bias score" (what I'm calling "asymmetry index"), as outlined in Graham et al., 2023. 
%   Input comes from ------------   [X, C, neuronList, numSess, cfg_out] = AHV_pearson_correlation_CONTRA_IPSI(cfg_in, tfilelist)   ---------------------

CWs = C.slopeCW;
CCWs = C.slopeCCW;
Smax = NaN(1,length(tfilelist));

    %% Which side has the greater value      1 = CW. 0 = CCW
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




















