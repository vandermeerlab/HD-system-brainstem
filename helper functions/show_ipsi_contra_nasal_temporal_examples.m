% cd('C:\Jeff\U01\datatouse\M212\M212-2021-07-21')  % random example session
function show_ipsi_contra_nasal_temporal_examples(startframe)
% 2025-04-03. JJS. This function displays eye tracking camera image with computed pupil position superimposed. 
% Inputs:   startfame -- which frame to start the video on. Find a good segement (during AHV), and divide time by sample step (.02s) to get frame #.
SSN = HD_GetSSN;
counter = 0;
if exist('events_ts.mat') == 2
    load('events_ts.mat');
else
    events_ts = LoadEvents([]);
end
index = strfind(events_ts.label, 'Starting Recording');
%             temp = cellfun(@isempty, index, 'UniformOutput', false);
for iCell = 1:length(index)
    isone(iCell) = ~isempty(index{iCell});
end
indextouse = find(isone);
if length(indextouse) == 1
    starttime = events_ts.t{indextouse};
else
    error('could not find start time for this session')
end

%% Get timestamps from the .smi file
[~, b, c] = fileparts(FindFile('*VT1.smi'));
fn = strcat(b,c);
tvec_raw = read_smi(fn);
tvec = tvec_raw - starttime;

%% Get the pupil trace dervied from Facemap
f = FindFiles('*VT1_proc.mat');
load(f{1}, 'pupil', 'rois');
pupilH = pupil{1}.com(:,2);
pupilV = pupil{1}.com(:,1);

plot(tvec, pupilH)

axis([843.0245  849.8004   46.3839  127.1283]);
vidObj = VideoReader(strcat(SSN, '-VT1.mp4'));
% t1 = 843.0245; % somewhat arbitrary starting point

% [value, index] = min(abs(tvec - t1));
% startframe = index;
% endframe = index + 100;

% startframe = 42232;
endframe = startframe + 500;

%% for loop to walk through video frames and display the associated pupil position in pixel coordinates

close all
% Create a figure window in full screen mode
hFig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

% Create axes to cover the entire figure
axes('Position', [0 0 1 1]);  % Full screen axes (normalized units)

% Loop through each frame of the video
for iframe = startframe:endframe
    counter = counter + 1;
    % Read the next frame
    frame = read(vidObj, iframe);
    %     frame = readFrame(vidObj);
    
    % Display the frame using imshow with specified axes
    imshow(frame, 'Parent', gca);  % Display on the current axes (gca)
    hold on
    
    first_fb = imlocalbrighten(frame);    % Brighten the image. Sometimes the image is dark and brightening it helps to pick out the new axis.
    clf; imagesc(first_fb); hold on
    
    disp(num2str(iframe))
    title(strcat(SSN, '--','Frame #', num2str(iframe)));
    xlabel(strcat('timestamp = ', num2str(tvec(iframe))))
    set(gca, 'FontSize', 20)
    
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
    %     plot(xNew, yNew, '.', 'Color', [.9 .9 .9], 'MarkerSize', .05) % plot all of the pupil positions for the session
    scatter(xNew, yNew, 'MarkerFaceAlpha', .01, 'MarkerEdgeAlpha', .01)
    
    % Plot pupil center of mass for this frame only
    %     plot(pupil{1}.com(iframe,2) + xmin, pupil{1}.com(iframe,1) + ymin, 'MarkerSize', 50, 'Color', 'r');  % WHY ISN'T plot working?!
    xval(counter) = pupil{1}.com(iframe,2) + xmin; % horizontal position of the pupil in pupil ROI frame of reference
    yval(counter) = pupil{1}.com(iframe,1) + ymin;
    scatter(xval(counter), yval(counter), 'LineWidth', 20, 'MarkerFaceColor', [1 1 0]) % not sure why this dot isn't yellow, but whatever
    text(400, 300, strcat('x-position (pixels) =', num2str(round(xval(counter),1))), 'FontSize', 35, 'Color', 'r')
    if counter ~= 1
        text(400, 350, strcat('delta x (pixels) =', num2str(round(xval(counter) - xval(counter - 1),1))), 'FontSize', 35, 'Color', 'r')
    end
    pause
    
    
end
















