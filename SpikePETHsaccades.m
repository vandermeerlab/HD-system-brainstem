function SpikePETHsaccades(S, saccade_input, varargin)

% 2021. JJS. 

process_varargin(varargin);

s = inputname(2);

if strcmp(s, 'nasalSaccades')
    titletouse = 'nasalSaccades';
elseif strcmp(s, 'temporalSaccades')
    titletouse = 'temporalSaccades';
else
    titletouse = [];
    warning('input not expected')
end

figure
[outputS, outputT, outputGau, pre_stim_mean, post_stim_mean] = SpikePETHvdm([], S, saccade_input);
subplot(2,1,1)
set(gca, 'FontSize', 20)
set(gca, 'XTick', [-0.5 -0.4 -0.3 -0.2 -0.1 0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1])
xlabel('')
title(S.label{1})

subplot(2,1,2)
set(gca, 'FontSize', 20)
set(gca, 'XTick', [-0.5 -0.4 -0.3 -0.2 -0.1 0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1])
xlabel('Time (sec)')
ylabel('Firing Rate', 'FontSize', 20)
title(titletouse)
c = axis;
axis([-0.5 1 c(3) c(4)])
