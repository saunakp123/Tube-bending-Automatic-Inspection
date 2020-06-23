% Chris Woodruff 6/25/19
% This program calibrates the camera based on a series of 15-20 images of
% the checkerboard calibration pattern. Read more at this link:
% https://www.mathworks.com/help/vision/ug/single-camera-calibrator-app.html
clc, clear all, close all;
% change to folder location of calibration images
images = imageDatastore('../Camera Calibration/');
% start calibration sequence
[imagePoints,boardSize] = detectCheckerboardPoints(images.Files(1:end));
% update squareSize if new pattern is printed
squareSize = 6.88; % millimeters
worldPoints = generateCheckerboardPoints(boardSize,squareSize);
% read each image in folder
I = imread(images.Files{end}); 
imageSize = [size(I,1) size(I,2)];
% identify Camera parameters
[cameraParams, ~, estimationErrors] = estimateCameraParameters(imagePoints,worldPoints,...
    'ImageSize',imageSize);
% display calibration results
figure; 
showExtrinsics(cameraParams, 'CameraCentric');
figure; 
showExtrinsics(cameraParams, 'PatternCentric');
figure; 
showReprojectionErrors(cameraParams);
displayErrors(estimationErrors, cameraParams);
% save camera parameters for image undistortion
save('../Camera Calibration/cameraParams.mat')
