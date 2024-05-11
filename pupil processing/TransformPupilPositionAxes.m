%
data = [1 1];

new_axis_start = [-2 0];
new_axis_end = [2 1];

%
subplot(221)
plot([new_axis_start(1) new_axis_end(1)], [new_axis_start(2) new_axis_end(2)], 'b', 'LineWidth', 2);
hold on;
plot(data(1), data(2), '.r', 'MarkerSize', 10)

%% translate to origin
data_r = data - new_axis_start;
new_axis_start_r = new_axis_start - new_axis_start;
new_axis_end_r = new_axis_end - new_axis_start;

%
subplot(222)
plot([new_axis_start_r(1) new_axis_end_r(1)], [new_axis_start_r(2) new_axis_end_r(2)], 'b', 'LineWidth', 2);
hold on;
plot(data_r(1), data_r(2), '.r', 'MarkerSize', 10)

%% projection
pr = (dot(data_r,new_axis_end_r)/norm(new_axis_end_r).^2)*new_axis_end_r;
plot(pr(1), pr(2), '.g', 'MarkerSize', 10)
axis equal;

%% move projected point back to original coordinates
pr_orig = pr + new_axis_start;

subplot(221)
plot(pr_orig(1), pr_orig(2), '.g', 'MarkerSize', 10)
axis equal;