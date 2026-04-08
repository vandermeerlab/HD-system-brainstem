% plot_session_FR_EV.m
% Full-session FR + Eye Velocity trace with stationary and laser shading.
% cd to session directory first, then run.

% --- Config ---
neuron_name = 'M284-2022-05-12-1-TT07_1';

% --- Load spikes ---
SSN = HD_GetSSN;
[S_all, cell_names] = LoadSpikesQuick();
idx = find_cell_index(cell_names, neuron_name, SSN);
myCell = SelectTS([], S_all, idx);

% --- Compute FR via MakeQfromS ---
cfg_Q = []; cfg_Q.dt = 0.001; cfg_Q.gausswin_sd = 0.05; cfg_Q.smooth = 'gauss';
Q = MakeQfromS(cfg_Q, myCell);

% --- Load eye velocity ---
sacc = load(strcat(SSN, '-saccades-edited.mat'), 'diffHdeg', 'diffH');
if isfield(sacc, 'diffHdeg'); EV = sacc.diffHdeg; else; EV = sacc.diffH; end

% --- Load stationary times ---
stat = load(strcat(SSN, '-AHV_StationaryTimes.mat'), 'STtstart', 'STtend');
st_starts = stat.STtstart(:);
st_ends   = stat.STtend(:);

% --- Check for laser ---
hasLaser = false;
try
    [~, ~, laser_on, laser_off, arraysize] = SortBrainstemEventLabels3;
    if ~isnan(laser_on); hasLaser = true; end
catch
end

% --- Plot (same as version that made the good figure) ---
figure('Position', [50 50 1800 500]); hold on
yyaxis left
plot(Q.tvec, Q.data./cfg_Q.dt, 'LineWidth', 1);
ylabel('FR (Hz)', 'FontSize', 18)
set(gca, 'FontSize', 14)

yyaxis right
plot(EV.tvec, EV.data, 'LineWidth', 0.5);
ylabel('Eye Vel (deg/s)', 'FontSize', 18)

xlabel('Time (s)', 'FontSize', 18)
title(neuron_name, 'FontSize', 20, 'Interpreter', 'none')
box off

% --- Now add shading AFTER axes are scaled ---
yyaxis left
c = axis;

% Stationary shading (light grey, behind data)
for iST = 1:length(st_starts)
    t0 = st_starts(iST);
    t1 = st_ends(iST);
    fill([t0 t1 t1 t0], [c(3) c(3) c(4) c(4)], ...
        [0.9 0.9 0.9], 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Laser shading (cyan)
if hasLaser
    for iL = 1:arraysize
        rectangle('Position', [laser_on(iL), c(3), laser_off(iL)-laser_on(iL), c(4)-c(3)], ...
            'FaceColor', [0 1 1], 'EdgeColor', [0 1 1]);
    end
end

% Replot data on top of shading
yyaxis left
plot(Q.tvec, Q.data./cfg_Q.dt, 'LineWidth', 1);
yyaxis right
plot(EV.tvec, EV.data, 'LineWidth', 0.5);
