function [optoPETHall, cfg] = plot_average_optotag_laser_PETH(optosessionstouse)

% inputs
%       myOptoCell = which ordered number a neuron is within its own session folder. For ex., ...TT01_1 would be a [1]. A later neuron would have a different, larger number.

myOptoCell = [1          2          1          1          1          4          2]; % 1st 7 opto cells
          
cfg_laser.dt = 0.1;
cfg_laser.doPlot = 0;
cfg_laser.doRaster = 0;
cfg_laser.doBar = 0;

if isempty(optosessionstouse)
    optosessionstouse = vertcat({'C:\Jeff\U01\datatouse\M271\M271-2021-08-30\M271-2021-08-30-TT02_1.t'    }, ...
        {'C:\Jeff\U01\datatouse\M282\M282-2022-02-10\M282-2022-02-10-TT03_2.t'    }, ...
        {'C:\Jeff\U01\datatouse\M282\M282-2022-02-15\M282-2022-02-15-TT07_1.t'    }, ...
        {'C:\Jeff\U01\datatouse\M391\M391-2022-12-07-1\M391-2022-12-07-1-TT01_1.t'}, ...
        {'C:\Jeff\U01\datatouse\M391\M391-2022-12-07-2\M391-2022-12-07-2-TT02_1.t'}, ...
        {'C:\Jeff\U01\datatouse\M402\M402-2022-11-02-2\M402-2022-11-02-2-TT02_3.t'}, ...
        {'C:\Jeff\U01\datatouse\M444\M444-2023-08-18\M444-2023-08-18-TT07_1.t'    });
end

doPlotIndividual = 0;
if doPlotIndividual
    cfg_laser.window = [-1.1 2]; %#ok<UNRCH>
    cfg_laser.dt = 0.01;
    cfg_laser.binsize = cfg_laser.dt; % used for gaussian kernal.  select a small bin size for good time resolution
    cfg_laser.doPlot = 1;
    cfg_laser.doRaster = 1;
    cfg_laser.doBar = 0;
end

doPlotAverage = 0;
FontSize = 12;

sessCounter = 0;
for iNeuron = 1:length(optosessionstouse)
    [path, neuron_to_use, ~] = fileparts(optosessionstouse{iNeuron});
    if strcmp(pwd, path) == 0
        pushdir(path);
        sessCounter = sessCounter + 1;
        EvalKeys
       
        sd = LoadSessionData([]); 
        myCell = SelectTS([], sd.S, myOptoCell(iNeuron));
        [~, cellID, ~] = fileparts(optosessionstouse{iNeuron})
        
        if doPlotIndividual
            [start_time, stop_time, laser_on, laser_off, arraysize, difflaser] = SortBrainstemEventLabels3; 
            dur = mode(laser_off(1:arraysize) - laser_on(1:arraysize));
            clf
            [~, ~, ~, ~, ~] = SpikePETH_either(cfg_laser, myCell, laser_on); hold on
            c = axis;
            rectangle(Position = [0, 0, dur, c(4)], FaceColor=[0 1 1], EdgeColor=[0 1 1])    % change the 3rd entry in rectangle to the diff of laser_off and laser_on
            [outputS, outputT, outputGau, outputIT, cfg] = SpikePETH_either(cfg_laser, myCell, laser_on);  % this is a hack to get the background color 'in back'. otherwise it occludes the spikes.
            line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'w')
            line([1 1], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'w')
            set(gca, 'YTick', [])
            set(gca, 'FontSize', FontSize)
            title('Laser PETH', 'FontSize', FontSize)
            set(gca, 'XTick', [-1 0 1 2])
            hold on
            xlabel('Time (s)')
            % add line for Ave Firing Rate, esp. when rasters are so dense that it turns black
            [outputS, outputT, outputGau, outputIT, cfg] = SpikePETH_either(cfg_laser, myCell, laser_on);
            m = histc(outputS, outputIT);
            title(cellID)
            yyaxis right
            plot(outputIT,m/cfg.dt/length(laser_on), 'LineWidth', 6, 'Color', 'r');
            ax = gca;
            ax.YAxis(2).Color = 'r'; 
            ylabel('Firing Rate (Hz)')
            disp('press any key to continue')
            pause            
        else
            [start_time, stop_time, laser_on, laser_off, arraysize, difflaser] = SortBrainstemEventLabels3; %#ok<*UNRCH>
            dur = mode(laser_off(1:arraysize) - laser_on(1:arraysize));
            [outputS, outputT, outputGau, outputIT, cfg] = SpikePETH_either(cfg_laser, myCell, laser_on);
            temp = histc(outputS, outputIT);
            optoPETHall(iNeuron,:) = temp/cfg.dt/length(laser_on);
        end
    end
end

if doPlotAverage
    
    
    
    
    
    
    
    
    
    
end

cfg.xbin = cfg.window(1):cfg.dt:cfg.window(2);













