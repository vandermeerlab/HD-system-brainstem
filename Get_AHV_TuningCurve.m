function Get_AHV_TuningCurve(AHVtsd, S, varargin)
%2020-01-14. JJS. Plots a tuning curve of firing rate for all cells versus the angular head velocity (AHV).
%   Detailed explanation goes here

plotRawTuning = 0;
extract_varargin;

% Ftsd = tsd(F.range, F.data);
TC = TuningCurves([],S,AHVtsd);

binsToUse = TC.min:(TC.max-TC.min)/(TC.nBin-1):TC.max;

fc = FindFiles('*.t', 'CheckSubdirs', 0);
numCells = size(TC.H,1);
figure
for iCell = 1:numCells;
    subplot(1,numCells,iCell)
    plot(-binsToUse, TC.H(iCell,:)./TC.Occ', 'LineWidth', 5);    % binsToUse gets reversed because we need to flip things about the y-axis. The way the encoder is set up with the
    % arduino, clockwise rotations yield positive values and CCW negative. But the convention in the field is that CCW turns are pos and CW turns negative.
    %     pause
    [a, b, c] = fileparts(fc{iCell}); %#ok<NASGU,ASGLU>
    title(b, 'FontSize', 18)
    set(gca, 'FontSize', 18)
    ylabel('Firing Rate', 'FontSize', 18)
    xlabel('AHV', 'FontSize', 18)
    %     clf
end

if plotRawTuning == 1;
    F = FiringRateNoSD(S, 0.1);
    Ftsd = tsd(F.range, F.data);
    plot(-AHV.data, Ftsd.data(AHV.range), '.')
    title(b, 'FontSize', 18)
    set(gca, 'FontSize', 18)
    ylabel('Firing Rate', 'FontSize', 18)
    xlabel('AHV', 'FontSize', 18)
end




