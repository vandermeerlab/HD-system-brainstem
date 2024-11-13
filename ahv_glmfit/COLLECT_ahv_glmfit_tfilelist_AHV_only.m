function [cfg_master, out] = COLLECT_ahv_glmfit_tfilelist_AHV_only(tfilelist)
% JJS. 2024-10-02. This version works from a tfile list of .t files (path included), instead of pushing in and out of each session folder to calculate all neurons.
%                  This way we can calculate for only NPH neurons, for example.
% This version plots all points in the AHV and pupil tuning curves (in addition to the mean)
%% collect data from all sessions
%%
cfg_master = []; % overall params
cfg_master.doColor = 1;
cfg_master.FontSize = 10;
cfg_master.doPlot = 0;
cfg_master.WriteFig = 1;
cfg_master.dt = 0.005;
cfg_master.maxlag = 200; % bins for use in saccade PETH
cfg_master.debug = 0;
cfg_master.tc_binEdges = -150:10:150;
cfg_master.pupil_tc_binEdges = -80:5:80; % for pupilX TC
cfg_master.doCat = 0;
cfg_master.doSave = 0;
cfg_master.doPlot = 0;

out = [];
for iNeuron = 1:length(tfilelist)
    [out{iNeuron}] = SESSION_ahv_glmfit_tfilelist_AHV_only(cfg_master, tfilelist{iNeuron});
    
end

%% ERROR with concatenating. Look at matt's original version.
% ahv_gain_fun = @(x) x.pca_sacc_both.rsq - x.pca_sacc.rsq;
% sacc_gain_fun = @(x) x.pca_sacc_both.rsq - x.ahv.rsq;
% rsq_fun = @(x) x.pca_sacc_both.rsq;
%
% ahv_gain = arrayfun(ahv_gain_fun, out);
% sacc_gain = arrayfun(sacc_gain_fun, out);
% rsq_all = arrayfun(rsq_fun, out);

cd('C:\Jeff\U01\datatouse');
if cfg_master.doSave == 0;
    disp('saving')
    save('GLM.mat')
end
%% Concatenate
ahv_true = [];
ahv_abs = [];
ahv_both = [];

if cfg_master.doCat
    for iNeuron = 1:length(out)
        ahv_true.rsq(iNeuron) = out{1,iNeuron}.ahv_true.rsq;
        ahv_true.p(iNeuron) = out{1,iNeuron}.this_m.Coefficients{3,4};
        ahv_true.t(iNeuron) = out{1,iNeuron}.this_m.Coefficients{3,3};
        
        ahv_abs.rsq(iNeuron) = out{1,iNeuron}.ahv_absolute.rsq;
        ahv_abs.p(iNeuron) = out{1,iNeuron}.this_m.Coefficients{4,4};
        ahv_abs.t(iNeuron) = out{1,iNeuron}.this_m.Coefficients{4,3};
        
        ahv_both.rsq(iNeuron) = out{1,iNeuron}.ahv_both.rsq;
        %         ahv_both.p(iNeuron) = out{1,iNeuron}.this_m.Coefficients{5,4};
        %         ahv_both.t(iNeuron) = out{1,iNeuron}.this_m.Coefficients{5,3};
    end
end


% PLOT
clf; doSave = 1;
if doSave; cd('D:\Jeff\U01\analysis\AHV GLM\TCs with Rsq'); end
doPause = 0; 
doPlot = 1;

if doPlot
    for iNeuron = 1:length(out)
        plot(out{1,iNeuron}.ahvscatterY, out{1,iNeuron}.ahvscatterX, '.')
        xlabel('AHV')
        ylabel('FR')
        title(num2str(iNeuron))
%         set(gca, 'FontSize', 26)
        y = ylim;
        text(-190, .95*(y(2)), strcat('true=', num2str(round(ahv_true.rsq(iNeuron),2))), 'FontSize', 20);
        text(-190, .87*(y(2)), strcat('abs=', num2str(round(ahv_abs.rsq(iNeuron),2))), 'FontSize', 20);
        c = axis;
        axis([-200 200 0 c(4)]);
        if doPause
            pause
        end
        if doSave
            disp('saving fig')
            %             savefig(num2str(iNeuron))
            saveas(gcf,strcat(num2str(iNeuron), '.png'))
        end
    end
end







end

