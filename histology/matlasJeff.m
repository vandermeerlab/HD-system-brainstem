function [] = matlasJeff

%% start with vandermeerlab path
% addpath(genpath('C:\Users\mvdm\Dropbox\projects\Jeff\matlas\'));
datapath = 'D:\Jeff\U01\histology\test folder';
cd(datapath);

%% load all images (currently anything in wd that ends in .png)
global cfg_master imglib curr_fig_no

cfg_master.markersize = 20;
cfg_master.alpha = 0.6;

pushdir(datapath);
imglib.fc = FindFiles('*.jpg'); % note that we'll probably end up with different sets of images (e.g. with and without labels), these should be named so that FindFiles can separate them
popdir;

nF = length(imglib.fc); if nF == 0, error('No .png files found.');end

% assemble main data strucure
for iF = 1:nF
    [~, imglib.fn{iF}, ~] = fileparts(imglib.fc{iF});
    imglib.imdata.(imglib.fn{iF}) = imread(imglib.fc{iF});
end

%% load histology database
imglib.histdata = readtable('brainstem.xlsx');

%% load image calibration data file if it exists
if exist('cal.mat', 'file')
   temp = load('cal.mat'); 
   imglib.cal = temp.cal; % cal.fn.xoffset, cal.fn.xmult, cal.fn.yoffset, cal.fn.ymult, (cal.fn.md5 once image checksums added)
else
   imglib.cal = [];
end

%% check if calibration data exists for each image file (need to implement checksums for extra robustness in case images get overwritten)
if isempty(imglib.cal)
    cal_keys  = {};
else
    cal_keys = fieldnames(imglib.cal);
end

to_calibrate = setxor(imglib.fn, cal_keys);

%% if not, prompt input to calibrate image data (click on two requested points) 
% then store offsets and multipliers to convert between pixels to mm, and
% save
for iF = 1:length(to_calibrate)
   
    close all;
    imh = imshow(imglib.imdata.(to_calibrate{iF}), 'Border', 'tight');
    hold on;
    
    disp('Click at ML -3.0, DV -5.0 (left of midline; use RIGHT axis for DV)');
    c1 = ginput; % first element is x, second is y coordinate (in image pixels)
    
    disp('Click at ML +3.0, DV -1.0 (right of midline; use RIGHT axis for DV):');
    c2 = ginput;
    
    xdiff = c2(1) - c1(1);
    imglib.cal.(to_calibrate{iF}).xoffset = mean([c1(1) c2(1)]); % pixel corresponding to midline
    imglib.cal.(to_calibrate{iF}).xmult = xdiff ./ 6; % pixels per mm
    
    ydiff = c2(2) - c1(2);
    imglib.cal.(to_calibrate{iF}).yoffset = c1(2) + ydiff * 1.25;
    imglib.cal.(to_calibrate{iF}).ymult = ydiff ./ 4; % pixels per mm
    
    vline(imglib.cal.(to_calibrate{iF}).xoffset);
    hline(imglib.cal.(to_calibrate{iF}).yoffset);
    pause(1);
    
    imglib.cal.(to_calibrate{iF}).ap = input('AP coordinate: ');
    
end

cal = imglib.cal; save('cal.mat', 'cal');

%% display first image
curr_fig_no = 1;

figure('KeyPressFcn', @fig_callback)
imshow(imglib.imdata.(imglib.fn{curr_fig_no}), 'Border', 'tight');
hold on;

%% can now use figure callbacks to do stuff (press 'h' for help)