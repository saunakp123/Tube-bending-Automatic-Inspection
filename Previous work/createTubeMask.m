function [mask] = createTubeMask(image,theta_bend,d_tube,offset,r_bend,xCenter,yCenter,length1,length2)
%CreateTubeMask(image,theta,d_tube,offset,r_bend,xCenter,yCenter)
%  Generate a background mask for first step of image filtering.
%
%  Inputs
%      image: reference image for size
%      theta: input bend angle
%      d_tube: tube diameter
%      offset: mask offset from expected tube edge
%      r_bend: bend radius
%      xCenter: bend centerpoint
%      yCenter: bend centerpoint
%      length1: first segment length
%      length2: second segment length

%  Output:
%      mask: background mask image
%
% create zero image to match input image size
imSize = size(image);
imageSizeX = imSize(2);
imageSizeY = imSize(1);
myImage1 = zeros(imageSizeY, imageSizeX, 'uint8');
% create curved portion of tube
% define inner and outer mask edge
r_in = r_bend - (d_tube/2 + offset);
r_out = r_bend + (d_tube/2 + offset);
R = [r_in:1:r_out];
for i=1:length(R)
    radius=R(i);
    theta = linspace(pi, pi+theta_bend, round(4* pi * radius)); % Define angles
    % Get x and y vectors for each point along the circumference.
    x = radius * cos(theta) + xCenter;
    y = radius * sin(theta) + yCenter;
% Write those (x,y) into the image with gray level 255.
    for k = 1 : length(x)
        row = round(y(k));
        col = round(x(k));
        myImage1(row, col) = 255;
    end
end
% figure;
% imshow(myImage1);
% create first straight section of tube
myImage2 = zeros(imageSizeY, imageSizeX, 'uint8');
y_start = imageSizeY - length1;
x_start = xCenter - r_out;
x_end = xCenter - r_in; 
for i=y_start:1:imageSizeY
    myImage2(i, x_start:x_end) = 255;  
end
% figure;
% imshow(myImage2);
% create second straight section of tube
myImage3 = zeros(imageSizeY, imageSizeX, 'uint8');
y_start = imageSizeY - length2;
x_start = xCenter - r_out;
x_end = xCenter - r_in; 
for i=y_start:1:imageSizeY
    myImage3(i, x_start:x_end) = 255;  
end
myImage3 = imtranslate(myImage3,[0, -yCenter]);
myImage3 = rotateAround(myImage3,yCenter,xCenter,-theta_bend*180/pi);
% figure;
% imshow(myImage3);
mask = myImage1 + myImage2 + myImage3; 