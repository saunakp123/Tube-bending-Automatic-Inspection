% Chris Woodruff 6/25/19
% This program develops a background mask for filtering out noise and other
% objects from the image background, based on the expected bend angle and
% length of the tube. 
clc, clear all, close all;
% Variables to be defined based on actual image size, x/y location of the
% center of the bend, bend radius, tube length (all in pixels)
imageSizeX = 1200;
imageSizeY = 1200;
myImage = zeros(imageSizeY, imageSizeX, 'uint8');
xCenter = 500;
yCenter = 600;
r_bend = 400;
d_tube = 50; 
offset = 30;
length1 = 600;
length2 = 200;
theta_bend = 3*pi/4;
mask = createTubeMask(myImage,theta_bend,d_tube,offset,r_bend,xCenter,yCenter,length1,length2);
figure;
imshow(mask);
axis('on', 'image');
% Plot crosshairs in the overlay at the center
hold on;
plot(xCenter, yCenter, 'r+', 'MarkerSize', 100);