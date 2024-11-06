function [tc_out_AHV] = senzai_functions_to_use(mouseIDs, varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
HD_TCbinCenters = -175:10:175;

doPlot = 1;
doSmooth = 0;
doSave = 1;
doPause = 0;
process_varargin(varargin)

if isempty(mouseIDs)
    mouseIDs = horzcat({'77c'},{'79c'},{'83b'},{'85b'},{'96d'},{'98b'});
end

for iMouse = 1:6
    IDtoUse = mouseIDs{iMouse};
    [~,~, ~, ~, HD_TCs{iMouse},~,~] = senzai_extract_spike_train(IDtoUse);
    [HD_neuronsToUse{iMouse}] = senzai_choose_HD_cells(HD_TCs{iMouse}, 10, 'doPlot', 0);
    [~, AHVtsd] = senzai_get_HD_AHV(IDtoUse);
    [tc_out_AHV{iMouse}] = senzai_getAHV_TCs(HD_neuronsToUse{iMouse}, S, AHVtsd, 'doPlot', 0);
end

%% Concatenate the mice AHV TCs
clear ahv_TC_all
binCenters = tc_out_AHV{1,1}.binCenters;
neuronCounter_overall = 0;
HD_TCs_all = [];
HD_zscore__all = [];
HD_pval__all = [];
for iMouse = 1:length(tc_out_AHV)
    for iNeuron = 1:size(tc_out_AHV{1,iMouse}.tc,1)
        neuronCounter_overall = neuronCounter_overall + 1;
        ahv_TC_all(neuronCounter_overall,:) = tc_out_AHV{1,iMouse}.tc(iNeuron,:);
        if sum(HD_neuronsToUse{iMouse} == iNeuron) == 1
            isHD(neuronCounter_overall) = 1;
        else
            isHD(neuronCounter_overall) = 0;
        end
    end
    HD_TCs_all = horzcat(HD_TCs_all, HD_TCs{iMouse}.Maps');
    HD_zscore__all = horzcat(HD_zscore__all, HD_TCs{iMouse}.z_ppln');
    HD_pval__all = horzcat(HD_pval__all, HD_TCs{iMouse}.pval_ppln');
end

%% For plotting ALL ahv tuning curves (n = 6 mice)
if doPlot
    for iNeuron = 1:size(ahv_TC_all,1)
        clf;
        f = figure(1);
        f.Position = [300 200 1080 800];
        %         if ismember(iNeuron, isHD)
        %             colortouse = 'r';
        %         else
        %             colortouse = 'b';
        %         end
        yyaxis left
        if doSmooth
            plot(binCenters, smoothdata(ahv_TC_all(iNeuron,:)), 'Color', 'b', 'LineWidth', 8)
        else
            plot(binCenters, ahv_TC_all(iNeuron,:), '.', 'Color', 'b', 'MarkerSize', 20)
        end
        title(num2str(iNeuron))
        c = axis;
        axis([c(1) c(2) 0 c(4)]);
        xlabel('AHV (deg/s)')
        ylabel('FR for AHV')
        text(-150, 2, strcat('HD zscore =', num2str(round(HD_zscore__all(iNeuron),1))), 'FontSize', 15);
        set(gca, 'FontSize', 26)
        
        yyaxis right
        plot(HD_TCbinCenters, smoothdata(HD_TCs_all{iNeuron}.rate))
        ylabel('FR for HD')
        c = axis;
        axis([c(1) c(2) 0 c(4)]);
        if doPause
            pause
        end
        if doSave
            disp('saving fig')
            cd('D:\Jeff\U01\analysis\senzai dataset\tuning curves dot png')
%             savefig(num2str(iNeuron))
            saveas(gcf,strcat(num2str(iNeuron), '.png'))
        end
    end
end




