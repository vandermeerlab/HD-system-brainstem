function [AHVtsd, orientationtouse] = get_AHVtsd(cfg_in)
%2023-07-19. JJS. 
%   Calculates a tsd of AHV values (degrees per second) for an individual session. 
%
% Inputs: 
%           cfg_in:     variable inputs
%           sd:         structure with session data, include spike trains [sd.S]
%Outputs:
%           AHVtsd:     tsd of AHV values (deg./s) from the platform encoder on the recording rig. 


cfg_def.subsample_factor = 10;      % how much to subsample the AHV data
cfg = ProcessConfig(cfg_def, cfg_in);

[csc_tsd, hd_tsd, samplingrate, dt] = GetOrientationValues([]); %#ok<ASGLU>
orientationtouse = downsampleOrientationValues(hd_tsd, cfg.subsample_factor);

window = 0.1;
postsmoothing = .05;

tic 
AHV = dxdt(orientationtouse.tvec, orientationtouse.data, 'window', window, 'postsmoothing', postsmoothing);
toc
AHV = -AHV; % THIS STEP IS NECESSARY BECAUSE dxdt GIVES VALUES THAT ARE CORRECT, BUT WITH A SIGN FLIP.
AHVtsd = tsd(orientationtouse.tvec, AHV);


end



% cfg_AHV = [];
% cfg_AHV.subsample_factor = 10;
% [AHV_tsd, tc_out] = AHV_tuning(cfg_AHV, S);
% AHV_dt = median(diff(AHV_tsd.tvec));