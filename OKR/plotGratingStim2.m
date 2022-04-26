function [] = plotGratingStim2(numReps, varargin)
%2022-03-09. JJS. 
%   Plots sinusoidal grating stimuli in alternating bouts of LEFTward and RIGHTward motion. Uses the Visual Stimulus Toolbox 
%   https://github.com/UCI-CARL/VisualStimulusToolbox 
connectToCheetah

num

process_varargin(varargin)

objLeft = GratingStim([200 400], 'k', 50, 180, [.1 .3], 1, 0);
objRight = GratingStim([200 400], 'k', 50, 0, [.1 .3], 1, 0);
for iRep = 1: numReps

fullfig;
h = axes; set(h, 'position', [0 0 1 1]);
NlxSendCommand('-PostEvent "LEFT grating" 128 000')
plot(objLeft);

fullfig;
h = axes; set(h, 'position', [0 0 1 1]);
NlxSendCommand('-PostEvent "RIGHT grating" 128 111')
plot(objRight);

fullfig;
h = axes; set(h, 'position', [0 0 1 1]);
NlxSendCommand('-PostEvent "LEFT grating" 128 000')
plot(objLeft);

fullfig;
h = axes; set(h, 'position', [0 0 1 1]);
NlxSendCommand('-PostEvent "RIGHT grating" 128 111')
plot(objRight);
end


% pause(0.00001);
% frame_h = get(handle(gcf),'JavaFrame');
% set(frame_h,'Maximized',1);
% set(gcf,'units','normalized','outerposition',[0 0 1 1])