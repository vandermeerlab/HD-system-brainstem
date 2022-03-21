function [] = plotGratingStim()
%2022-03-09. JJS. 
%   Plots sinusoidal grating stimuli in alternating bouts of LEFTward and RIGHTward motion. Uses the Visual Stimulus Toolbox 
%   https://github.com/UCI-CARL/VisualStimulusToolbox 
connectToCheetah

% import java.awt.Robot;
% import java.awt.event.*;
% mouse = Robot;
% mouse.mouseMove(100, 100);
% mouse.mousePress(InputEvent.BUTTON2_MASK);    %left click press
% mouse.mouseRelease(InputEvent.BUTTON2_MASK);  %left click release

objLeft = GratingStim([200 400], 'k', 50, 180, [.1 .3], 1, 0);
objRight = GratingStim([200 400], 'k', 50, 0, [.1 .3], 1, 0);

fullfig;
h = axes; set(h, 'position', [0 0 1 1]);
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
NlxSendCommand('-PostEvent "LEFT grating" 128 000')
plot(objLeft);

import java.awt.Robot;
import java.awt.event.*;
mouse = Robot;
mouse.mouseMove(800, 800); % relates to my personal location of open program instance
mouse.mousePress(InputEvent.BUTTON1_MASK); % actual left click press
pause(0.001);
mouse.mouseRelease(InputEvent.BUTTON1_MASK); % actual left click release

fullfig;
h = axes; set(h, 'position', [0 0 1 1]);
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
NlxSendCommand('-PostEvent "RIGHT grating" 128 111')
plot(objRight);

fullfig;
h = axes; set(h, 'position', [0 0 1 1]);
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
NlxSendCommand('-PostEvent "LEFT grating" 128 000')
plot(objLeft);

fullfig;
h = axes; set(h, 'position', [0 0 1 1]);
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
NlxSendCommand('-PostEvent "RIGHT grating" 128 111')
plot(objRight);
end
