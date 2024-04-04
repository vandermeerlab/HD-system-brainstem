nSpikes = 1000; 
t_start = 0; t_end = 100;

fake_S = ts;
fake_S.t{1} = sort(t_start + (t_end - t_start) .* rand(nSpikes, 1));

fprintf('Expected mean FR: %.2f Hz\n', nSpikes ./ (t_end - t_start));
fprintf('Observed mean FR: %.2f Hz\n', length(fake_S.t{1}) ./ (t_end - t_start));

% tuning variable 1 (random uniform)
Fs = 50;
tvec = t_start:1/Fs:t_end;
fake_tsd = tsd(tvec, rand(size(tvec)));

cfg_tc = [];
cfg_tc.occ_dt = 1/Fs;
tc = TuningCurves(cfg_tc, fake_S, fake_tsd);

fprintf('Total occupancy %.2f s, expected %.2f +/- %.2f (1 sample)\n', sum(tc.occ_hist), (t_end - t_start), 1./Fs);

tc_fr = sum(tc.occ_hist .* tc.tc) ./ sum(tc.occ_hist); % average firing rate, weighted by occupancy
fprintf('Weighted mean FR in TC: %.2f Hz\n', tc_fr); % should equal true firing rate

% now cut out some data
remove_start = 5:10:95; remove_end = remove_start + 1; % intervals to remove

[~, keep] = restrict(fake_tsd, remove_start, remove_end);
fake_tsdR = fake_tsd; fake_tsdR.tvec = fake_tsdR.tvec(~keep); fake_tsdR.data = fake_tsdR.data(~keep);

[~, keep] = restrict(fake_S, remove_start, remove_end);
fake_SR = fake_S; fake_SR.t{1}(keep{1}) = [];
fprintf('Observed mean FR after restrict: %.2f Hz\n', length(fake_SR.t{1}) ./ ((t_end - t_start) - sum(remove_end - remove_start)));

tc = TuningCurves(cfg_tc, fake_SR, fake_tsdR);
fprintf('Total occupancy %.2f s, expected %.2f +/- %.2f (1 sample)\n', sum(tc.occ_hist), ((t_end - t_start) - sum(remove_end - remove_start)), 1./Fs);

tc_fr = sum(tc.occ_hist .* tc.tc) ./ sum(tc.occ_hist); % average firing rate, weighted by occupancy
fprintf('Weighted mean FR in TC: %.2f Hz\n', tc_fr); % should equal true firing rate
