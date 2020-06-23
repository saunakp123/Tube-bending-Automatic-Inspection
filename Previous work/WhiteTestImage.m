% Chris Woodruff 6/25/19
% This program runs the complete algorithm for calculating the bend angle
% by binarizing the image, identifying the centerline, and fitting a line
% to the centerline, taking the atan of the slope as the bend angle. 
%% Load Images
clc, clear all, close all;
% load image of tube bend
im = imread('./Images 53119/Part 2/Bend 2/white.png');
% load 'zero' image with no tube for background comparison
imz = imread('./Images 53119/whitez.png');
% load saved camera parameters
load './Camera Calibration/cameraParams.mat'
% display raw image
figure;
imshow(im);
%% Image Filtering
% undistort image and 'zero' image using calibrated camera parameters
[white_u,newOrigin] = undistortImage(im,cameraParams);
[white_uz,newOrigin] = undistortImage(imz,cameraParams);
% subtract 'zero' image from original
im=white_uz-white_u;
figure;
imshow(im);
% binarize image
im2=imbinarize(im);
figure; imshow(im2);
% fill in any holes caused by noisy image
bw2=imfill(im2,'holes');
figure;imshow(bw2);
[n_x n_y] = size(bw2);
xbar = 0;
bw3=0*bw2;
k=0;
% identify centerline by slicing each row of pixels
for i=n_x:-1:1
    if(sum(bw2(i,:))>180)
        k = k+1;
        xbar = (bw2(i,:)*(1:n_y)'/sum(bw2(i,:)));
        CL(k,:)= [i xbar];
        bw3(i,(round(xbar)-2):(round(xbar)+2)) = 1;
    end
end

% use this code if slicing the image by columns instead of rows
% for i=1:1650
%     if(sum(bw2(:,i))>180)
%         k = k+1;
%         SUM(k)=sum(bw2(:,i));
%         ybar = ((1:n_x)*bw2(:,i)/sum(bw2(:,i)));
%         N(k) = sum(bw2(:,i));
%         CL(k,:)= [ybar i];
%         bw3(round(ybar),i) = 1;
%     end
% end
CL = sortrows(CL,1);
% display image centerline
figure;
imshow(~bw3);
figure;
imshow(~bw3.*bw2);

%% derivative analysis
X = CL(:,2);
Y = CL(:,1);
figure;
plot(Y,X);
% hold on;
% vline([580 1280]);
% % subplot(211);

% fit spline to centerline data to smooth out pixel noise
Y2 = linspace(Y(1),Y(end),10000);
[fitresult, gof] = TubeSplineSmoothing(Y, X);
X2 = feval(fitresult,Y2);
% save tubeimagedata.mat

% calculate forward difference slope to identify tangent points
a = length(Y2);
for i = 1:(length(Y2)-1)
    dx(i) = (X2(i+1)-X2(i))/(Y2(i+1)-Y2(i));
end
figure;
plot(Y2(1:a-1),dx);
xlabel('y');
ylabel('dx/dy');
title('forward difference');

% manually input y coordinates of tangent break points
prompt = 'Input Break points [p1 p2 p3 p4 p5 p6] ';
a = input(prompt) %[300 600 700 1050 1450 2000]

% fit lines to spline points and calculate slope
for i=1:length(a)
    [d, ix] = min(abs(Y2-a(i)));
    aind(i)=ix;
end
a=aind;
[fitresult1, gof1] = linearTubeFit(Y2(a(1):a(2)), X2(a(1):a(2)));
[fitresult2, gof2] = createTubeFit(Y2(a(3):a(4)), X2(a(3):a(4)));
[fitresult3, gof3] = linearTubeFit(Y2(a(5):a(6)), X2(a(5):a(6)));
% pp = ppfit(sqrt(Y),X,3,[3 3 3]);

figure;
xlim([0 n_x]);
ylim([0 n_y]);
xlabel('y');
ylabel('x');
hold on;
plot(fitresult1);
plot(fitresult2);
plot(fitresult3);
plot(Y,X,'x');
% ylim([800 2000]);

% compare results to measured values
theta1 = atand(fitresult1.a)
theta_true = 50.47
err = theta_true - theta1

theta2 = atand(fitresult3.a)
theta_true = 0;
err = theta_true - theta2