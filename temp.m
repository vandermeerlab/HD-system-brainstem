%% Choose neurons with minimum FR and Strong HD tuning   [M77_neuronsToUse]
% Find the PFD (peak) for each neuron that has significant HD tuning
HD_z_thresh = 50;  % z-score threshold. ***Not sure how Yuta calculates the z-score. Check on this.
% Peak firing rate criterion
for iC = 1:length(M77_TCs.Maps)
    M77_maxFR(iC) = max(M77_TCs.Maps{iC,1}.rate);
    M77_maxFRsmooth(iC) = max(smoothdata(M77_TCs.Maps{iC,1}.rate));
end
FRthresh = 10; % in Hz
M77_peakFRtoUse = M77_maxFR > FRthresh;
M77_peakFRtoUseIDs = find(M77_peakFRtoUse); M77_peakFRtoUseIDs = M77_peakFRtoUseIDs';
sum_M77_peakFRtoUse = sum(M77_peakFRtoUse); fraction_M77_peakFRtoUse = sum_M77_peakFRtoUse/length(M77_TCs.Maps);

M77_HDC_logical = M77_TCs.z_ppln > HD_z_thresh;
M77_HDC_IDs = find(M77_HDC_logical);
M77_sum_HDCs = sum(M77_HDC_logical);

M77_neuronsToUse = intersect(M77_peakFRtoUseIDs, M77_HDC_IDs);
sum_M77_HDC_minFR = length(M77_neuronsToUse);

% clf; 
hold on; set(gca, 'FontSize', 24)
doNorm = 0;
if doPlot
    for iC = 1 : sum_M77_HDC_minFR
        clf
        IDtoUse(iC) = M77_neuronsToUse(iC);
        disp(num2str(IDtoUse(iC)))
        normIC(iC,:) = (M77_TCs.Maps{iC,1}.rate - min(M77_TCs.Maps{iC,1}.rate)) ./ (max(M77_TCs.Maps{iC,1}.rate) - min(M77_TCs.Maps{iC,1}.rate));
        M77_TCs.Maps{iC,1}.rateSmoothed = smoothdata(M77_TCs.Maps{iC,1}.rate);
        if doNorm
            plot(TCbinCenters, smoothdata(normIC(iC,:)));
            ylabel('normalized FR (Hz)')
        else
            plot(TCbinCenters, smoothdata(M77_TCs.Maps{IDtoUse(iC),1}.rate));
            line([halfwidth(iC,1) halfwidth(iC,1)], [0 M77_TCs_smoothed(iC, index1(iC))], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
            line([halfwidth(iC,2) halfwidth(iC,2)], [0 M77_TCs_smoothed(iC, index2(iC))], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
            ylabel('FR (Hz)')
        end
        title(strcat('M77 cell num', num2str(IDtoUse)))
        xlabel('HD (deg)')
        pause
    end
end
% c = axis;
% axis([0 365 c(3) c(4)]);