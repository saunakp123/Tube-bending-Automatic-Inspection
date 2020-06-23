% Chris Woodruff 6/25/19
% This program runs the complete algorithm for calculating the bend angle
% by binarizing the image, identifying the centerline, and fitting a line
% to the centerline, taking the atan of the slope as the bend angle. 
%% Load Images
clc, clear all, close all;
% load image of tube bend
im = imread('D:\Personal\UW\BARC\Part1\Bend 2\white1.jpeg');
% load 'zero' image with no tube for background comparison
imz = imread('D:\Personal\UW\BARC\Part1\whitez.jpeg');
% load saved camera parameters
load 'D:\Personal\UW\BARC\Checker\cameraParams.mat'
% display raw image
% figure;
% imshow(im);
%% Image Filtering
% undistort image and 'zero' image using calibrated camera parameters
[white_u,newOrigin] = undistortImage(im,cameraParams);
[white_uz,newOrigin] = undistortImage(imz,cameraParams);
% subtract 'zero' image from original
im=white_uz-white_u;
% figure;
 imshow(im);
% binarize image
 im2=imbinarize(im);
% figure; imshow(im2);
% fill in any holes caused by noisy image
bw2=imfill(im2,'holes');
figure;imshow(bw2);
[n_x n_y] = size(bw2);
xbar = 0;
bw3=0*bw2;
k=0;
l=0;
leftedge = 0*bw2 ;
rightedge = leftedge ;
%% Edge detection
% identifying edges 
for i=1:n_x 
    for j = 1:n_y-1
        if bw2(i,j)==0 && bw2(i,j+1)==1
            k = k+1 ;
            leftedge(i,j+1) = bw2(i,j+1) ;
            leftedgearr(k,:) = [j+1,i] ;
        end
        if bw2(i,j)==1 && bw2(i,j+1)==0
            l = l+1 ;
            rightedge(i,j) = bw2(i,j) ;
            rightedgearr(l,:) = [j,i] ;
            continue ;
        end
    end
end
%imshow(rightedge) 
%figure() ;
%imshow(leftedge)
edge = leftedge + rightedge ;
figure()
imshow(edge)
%% Angle calculation
%Left edge angle 
Yl = leftedgearr(:,2);
Xl = leftedgearr(:,1);
Yl_1 = leftedgearr(1:31,2);
Xl_1 = leftedgearr(1:31,1);
Yl_2 = leftedgearr(1091:1121,2);
Xl_2 = leftedgearr(1091:1121,1);
Yr = rightedgearr(:,2);
Xr = rightedgearr(:,1);
Yr_1 = rightedgearr(200:220,2);
Xr_1 = rightedgearr(200:220,1);
Yr_2 = rightedgearr(1091:1111,2);
Xr_2 = rightedgearr(1091:1111,1);
figure;
plot(Yl_1,Xl_1,'LineWidth',1.5)
axis equal
hold on
plot(Yl_2,Xl_2,'LineWidth',1.5)
hold on ;
plot(Yl,Xl);
axis equal
hold on
plot(Yr,Xr) ;
legend('1','1_1','0','2') ;
%Curvefitting to get slope of left edge line(before smoothing)
p = polyfit(Yl_1,Xl_1,1) ;
p1 = polyfit(Yr_2,Xr_2,1) ; % using right edge horizontal points 
angle = atan(p(1))*180/pi ;
angle1 = atan(p1(1))*180/pi ;
left_angle_presmooth = 180+angle-angle1
actual = 125.23 
error1 = abs(left_angle_presmooth-actual) 
%Curvefitting to get slope of right edge line(before smoothing)
p = polyfit(Yr_1,Xr_1,1) ;
p1 = polyfit(Yr_2,Xr_2,1) ;
angle = atan(p(1))*180/pi ;
angle1 = atan(p1(1))*180/pi ;
right_angle_presmooth = 180+angle-angle1
actual 
error2 = abs(right_angle_presmooth-actual)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fit spline to left edge to smooth out pixel noise
  Yl1 = linspace(Yl(1),Yl(end),10000);
  Yl1 = Yl1' ;
  [fitresult, gof] = smoothing(Yl, Xl);
  Xl1 = feval(fitresult,Yl1);
  % % fit spline to right edge to smooth out pixel noise
  Yr1 = linspace(Yr(1),Yr(end),10000);
  Yr1 = Yr1' ;
  [fitresult, gof] = smoothing(Yr, Xr);
  Xr1 = feval(fitresult,Yr1);
  %Plotting after smoothing
plot(Yl1,Xl1);
axis equal
hold on
plot(Yr1,Xr1) ;
 % Obtaining few points from the smooth tubeline to calculate angle 
Yl_1 = Yl1(1:30);
Xl_1 = Xl1(1:30);
Yl_2 = Yl1(5350:5450);
Xl_2 = Xl1(5350:5450);
Yr_1 = Yr1(920:950);
Xr_1 = Xr1(920:950);
Yr_2 = Yr1(5350:5450);
Xr_2 = Xr1(5350:5450);
  %Curvefitting to get slope of left edge line(after smoothing)
 p = polyfit(Yl_1,Xl_1,1) ;
 p1 = polyfit(Yr_2,Xr_2,1) ;
 angle = atan(p(1))*180/pi ;
 angle1 = atan(p1(1))*180/pi ;
 left_angle_smooth = 180+angle-angle1
 actual 
 error3 = abs(left_angle_smooth-actual)
 %Curvefitting to get slope of right edge line(after smoothing)
 p = polyfit(Yr_1,Xr_1,1) ;
 p1 = polyfit(Yr_2,Xr_2,1) ;
 angle = atan(p(1))*180/pi ;
 angle1 = atan(p1(1))*180/pi ;
 right_angle_smooth = 180+angle-angle1
 actual 
 error4 = abs(right_angle_smooth-actual)

% % calculate forward difference slope to identify tangent points
a = length(Yl);
for i = 1:(length(Yl)-1)
    dx(i) = (Xl(i+1)-Xl(i))/(Yl(i+1)-Yl(i));
end
figure;
plot(Yl(1:a-1),dx);
xlabel('y');
ylabel('dx/dy');
title('forward difference');
% 
% % manually input y coordinates of tangent break points
% prompt = 'Input Break points [p1 p2 p3 p4 p5 p6] ';
% a = input(prompt) %[300 600 700 1050 1450 2000]
% 
% % fit lines to spline points and calculate slope
% for i=1:length(a)
%     [d, ix] = min(abs(Y1_1-a(i)));
%     aind(i)=ix;
% end
% a=aind;
% [fitresult1, gof1] = linearTubeFit(Y1_1(a(1):a(2)), X1_1(a(1):a(2)));
% [fitresult2, gof2] = createTubeFit(Y1_1(a(3):a(4)), X1_1(a(3):a(4)));
% [fitresult3, gof3] = linearTubeFit(Y1_1(a(5):a(6)), X1_1(a(5):a(6)));
% % pp = ppfit(sqrt(Y),X,3,[3 3 3]);
% 
% figure;
% xlim([0 n_x]);
% ylim([0 n_y]);
% xlabel('y');
% ylabel('x');
% hold on;
% plot(fitresult1);
% plot(fitresult2);
% plot(fitresult3);
% plot(Y,X,'x');
% % ylim([800 2000]);
% 
% compare results to measured values
% theta1 = atand(fitresult1.a)
% theta_true = 50.47
% err = theta_true - theta1
% 
% theta2 = atand(fitresult3.a)
% theta_true = 0;
% err = theta_true - theta2