function [xNew, yNew, prX, prY, xT, yT, Hnew, pupil, rois] = recalculatePupilPosition(LeftEye, RightEye, first_f, cfg_in)

% JJS. 2024-04-29.
% This function recalculates the pupil position for each video frame into new coordinates, based on the axis defined by the left and right edges
% of the eye.
%   Inputs:     The coordinate inputs are obtained by using a different function to load the video file ______, and then using ginput to manually define the vertices and co-vertices of the eye (approximated as an ellipse).
%               LeftEye     - Pixel coordinates for the left edge of the eye. Note that this does not refer to the Left Eye. This refers to the left edge of the eye.
%                             Leave this variable as [] in the command line if you don't have it or want to regenerate it.
%               first_f     - first frame of the video file. For the current headfixed setup in B99, this should be a 480x640x3 uint8 array.
%               pupil       - part of the output from Facemap, which is stored in the '_proc.mat' file. pupil is a 1x1 cell array with one key value, the center of mass ('com'), or position of the pupil. nsamples x 2 vector.
%
%               *Note*:       There are two frames of reference with regard to the camera data. Each is expressed in pixel units (not degrees of visual angle).
%                           a) The axes for the camera image of the mouse (wide field of view). The size is 640x480 pixels. Note that on the y-axis, the values are descending from top left to bottom left.
%
%                           b) The axes for the pupil region of interest (ROI), which is rectangular area that is parallel to the wide field of view, but is user-defined and restricted to a small, subsection (eye only; actually, just part of the eye)
%                                   pupil{1,1}.com(:,2) [x-values]; These are calculate by Facemap and stored in the strcat(SSN,'_proc.mat') file in each session folder.
%                                   pupil{1,1}.com(:,1) [y-values];
%   Outputs:
%               xNew        - x-coordinates for pupil position, in the wide field of view frame of reference
%               yNew        - y-coordinates for pupil position, in the wide field of view frame of reference

%               prX         - the x-position of the pupil, in the new (eye-centered) axis
%               prY         - the y-position of the pupil, in the new (eye-centered) axis

%               xT          - the x-position of the pupil, translated to the original (camera) axis
%               yT          - the y-position of the pupil, translated to the original (camera) axis

% 2024-05-01.   Changed so that the top and bottom eye coordinates are eliminated from the ginput step. Instead, the minor axis of the pupil ROI is defined as 90 degrees (orthogonal) to the major axis.
% 2024-05-09.   Incorporated Mvdm's code from 'TransformPupilAxes.m' (my name for the script) in order to do the projection to the new axes accurately. This strategy uses the dot product.
%               I have left 'my' version of the position calculation at the bottom of this function, but commented out (for reference). It produced a nearly identical, but flipped (about the x-axis) version of the original pupil trace.
% 2024-05-10.   Implemented Mvdm's script and the output looks accurate now. Eye position values are different from the original trace, but not by much unless the angle btwn old and new axes is steep.
%
%               To do: 
%                           (1) Plotting all of the new data points and lines is prohibitively slow. It only really works if you wait an extremely long 
%                           time or run through just a few examples in debug mode. 
%                           (2) Need to include an option to scroll ahead to find a good frame if the mouse is blinking in the first frame of video. 
%                           (3) add as an output the angle between the camera horizontal and the new eye axis
%
cfg_def.doPlot = 0;
cfg_def.MarkerSize = .05;
cfg_master = ProcessConfig2(cfg_def, cfg_in);

%% Load the pupil data
SSN = HD_GetSSN;
fname = strcat(SSN, '-VT1_proc.mat');
if exist(fname)
    load(fname, 'pupil', 'rois')
else
    error('could not find facemap pupil data')
end
%% Load the video
if ~isempty('first_f')
    vname = strcat(SSN, '-VT1.mp4'); % filename of the video
    disp('Loading video...')
    vidObj = VideoReader(vname); % Load it into Matlab.
    first_f = read(vidObj, 1); % Read first frame. Usualy the first frame should be usable. If the mouse is blinking, however, it won't be. ***ADD***Option to scroll forward.
end
first_fb = imlocalbrighten(first_f);    % Brighten the image. Sometimes the image is dark and brightening it helps to pick out the new axis.
clf; imagesc(first_fb); hold on % This is the wide field of view image upon which everything else will be plotted.
%% Manually select the eye's corners.
if ~isempty('LeftEye') && ~isempty('RightEye')  % if these values were part of the input, skip this step
    disp('select the eyes corners, GOING IN THIS ORDER: Left, Right, enter')
    g = ginput; LeftEye = [g(1,1) g(1,2)]; RightEye = [g(2,1) g(2,2)];
end
xSPAN = abs(RightEye(1) - LeftEye(1)); % This is the width of the eye, if we end up needing that later.
if cfg_master.doPlot == 1
    plot(LeftEye(1), LeftEye(2), 'c.', 'MarkerSize', 40) % plot the endpoints of the major axis
    plot(RightEye(1), RightEye(2), 'c.', 'MarkerSize', 40)
end
%% Draw the ROI on top of the eye image
for iROI = 1:length(rois)
    a(iROI) = strcmp(rois{1,iROI}.rtype, 'pupil');
end
found = find(a);    % this is the cell array in the variable 'rois' that corresponds to the pupil ROI.
xrange = double(rois{1,found}.xrange); % xrange is a consecutive list of integers from xmin to xman. convert from uint8 to double
yrange = double(rois{1,found}.yrange); % yrange is a consecutive list of integers from ymin to yman. convert from uint8 to double
xspan = xrange(end) - xrange(1);  % length of the ellipse (x-axis). The ellipse major axis is always parallel to the wider field of view reference frame.
yspan = yrange(end) - yrange (1); % height of the ellipse (y-axis). The ellipse is usually longer in the x-dimension than the y-dimension.
xmin = xrange(1); xmax = xrange(end);
ymin = yrange(1); ymax = yrange(end);
[~, majorIndex] = max([xspan yspan]); % Identify the major axis of the ellipse. The pupil RIO is an ellipse. We will also plot the rectangle that bounds the ellipse.
if majorIndex == 1
    a = xspan/2;
    b = yspan/2;
    vertex1 = [xmin (ymin + yspan/2)]; % Identify the vertices of the major axis of the ellipse.
    vertex2 = [xmax (ymin + yspan/2)];
elseif majorIndex == 2
    a = yspan/2;
    b = xspan/2;
    vertex1 = [(xmin + xspan/2) ymin]; % Identify the vertices of the major axis of the ellipse.
    vertex2 = [(xmin + xspan/2) ymax];
else
    error('could not determine ellipse parameters')
end
x2 = vertex2(1); x1 = vertex1(1); % rename
y2 = vertex2(2); y1 = vertex1(2);
%% Calculate values that enable plotting
pupilROI_x  = [xmin xmax]; % x-coordinates for the rectangle that is the pupil ROI
pupilROI_y  = [ymin ymax]; % y-coordinates for the rectangle that is the pupil ROI

pupilX = pupil{1,1}.com(:,2); % x-coordinates for pupil position, in the whole camera field of view frame of reference
pupilY = pupil{1,1}.com(:,1); % y-coordinates for pupil position, in the whole camera field of view frame of reference
% pupilXwidth = pupilX./xSPAN;  % This is the x-coordinate for pupil position, expressed in fraction of an eye width.

xNew = pupilX + xmin; % x-coordinates for pupil position, in the wide field of view frame of reference
yNew = pupilY + ymin; % y-coordinates for pupil position, in the wide field of view frame of reference% e = sqrt(1 - b^2/a^2);  % e is the eccentricity of the ellipse. Did not end up needing this.

t = linspace(0,2*pi);
X = a*cos(t); Y = b*sin(t);
w = atan2(y2-y1,x2-x1); x = (x1+x2)/2 + X*cos(w) - Y*sin(w); y = (y1+y2)/2 + X*sin(w) + Y*cos(w); % Calculate the ellipse. Source?
%% Plot the pupil ROI ellipse
plot(vertex1(1), vertex1(2), 'c.', 'MarkerSize', 40) % plot the endpoints of the major axis
plot(vertex2(1), vertex2(2), 'c.', 'MarkerSize', 40)
plot(x,y,'c-', 'LineWidth', 5) % plot the ellispe
drawrectangle('Position',[pupilROI_x(1), pupilROI_y(1), xspan, yspan], 'StripeColor','r'); % draw a rectangle indicating the bounds of the pupil ROI from Facemap
hold on
plot(xNew, yNew, 'r.', 'MarkerSize', .05) % plot all of the pupil positions for the session
line([LeftEye(1) RightEye(1)], [LeftEye(2) RightEye(2)], 'Color', 'y', 'LineWidth', 5)
set(gca, 'FontSize', 50)
c = axis;
axis([c(1) c(2) -50 c(4)]); hold on
disp('Make sure everything looks good. Press enter to proceed.'); pause;
%% This begins Mvdm code section (adapted)
for iframe = 1 : length(pupil{1,1}.com)
    disp(num2str(iframe))
    data = [xNew(iframe) yNew(iframe)];
    new_axis_start = [LeftEye(1) LeftEye(2)];
    new_axis_end = [RightEye(1) RightEye(2)];
    %% Translate to new origin
    data_r = data - new_axis_start;
    new_axis_start_r = new_axis_start - new_axis_start;
    new_axis_end_r = new_axis_end - new_axis_start;
    %% projection
    pr = (dot(data_r,new_axis_end_r)/norm(new_axis_end_r).^2)*new_axis_end_r;
    % move projected point back to original coordinates
    pr_orig = pr + new_axis_start;
    %% Rename
    xT(iframe) = pr_orig(1); % xT stands for the x-position of the pupil, translated to the wider camera axis
    yT(iframe) = pr_orig(2); % xT stands for the y-position of the pupil, translated to the wider camera axis
    prX(iframe) = pr(1); % this is the x-position of the re-centered pupil point, in the original coordinate frame
    prY(iframe) = pr(2); % this is the y-position of the re-centered pupil point, in the original coordinate frame
    Hnew(iframe) = sqrt((prX(iframe) - 0)^2 + (prY(iframe) - 0)^2); % THIS IS THE OUTPUT OF INTEREST. It is the magnitude of the line that defines the horizontal pos. of the pupil one the new axis.
    
    if cfg_master.doPlot == 1
        %% Plot the new eye axis on a separate plot
        plot([new_axis_start(1) new_axis_end(1)], [new_axis_start(2) new_axis_end(2)], 'c', 'LineWidth', 4); hold on
        plot(data(1), data(2), '.r', 'MarkerSize', cfg_master.MarkerSize)
        plot([new_axis_start_r(1) new_axis_end_r(1)], [new_axis_start_r(2) new_axis_end_r(2)], 'c', 'LineWidth', 4);
        plot(data_r(1), data_r(2), '.r', 'MarkerSize', cfg_master.MarkerSize)
        line([0 pr(1)], [0 pr(2)], 'Color', 'y', 'LineWidth', 5)  % this line should start at origin and extend so as to meet the line coming down from the new pupil position at right angles
        line([0 Hnew(iframe)], [0 0], 'Color', 'g', 'LineWidth', 5)
    end
    clear pr
end



%% This section was from my original version where I used trig to try to get the new values. Did not work. 
% for iframe = 1 : length(pupil{1,1}.com)
%     xP = xNew(iframe);
%     yP = yNew(iframe);
%     d(iframe) = sqrt((xP - origin_x)^2 + (yP - origin_y)^2);
%     adj(iframe,1) = cos(90) * d(iframe);
%     opp(iframe,1) = sin(90) * d(iframe);
% end
%% misc
% [origin_x, origin_y] = linexline([LeftEye(1) RightEye(1)], [LeftEye(2) RightEye(2)], [TopEye(1) BottomEye(1)], [TopEye(2) BottomEye(2)] , 1); % this function finds the intersection point of two lines (our new origin)
%         axis([pupilROI_x(1) (pupilROI_x(1) + xspan) pupilROI_y(1) (pupilROI_y(1) + yspan)]) % restrict the image just to the pupil ROI
% axis([pupilROI_x(1) (pupilROI_x(1) + xspan) pupilROI_y(1) (pupilROI_y(1) + yspan)]) % restrict the image just to the pupil ROI

