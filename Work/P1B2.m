% Chris Woodruff 6/25/19
% This program runs the complete algorithm for calculating the bend angle
% by binarizing the image, identifying the centerline, and fitting a line
% to the centerline, taking the atan of the slope as the bend angle. 
% P1B1 107.26
% P1B2 125.23
%% Load Images
clc, clear all, close all;
% load image of tube bend
im = imread('whiteP1B2.png');
% load 'zero' image with no tube for background comparison
imz = imread('whitez.png');
% load saved camera parameters
load './Previous work/Camera Calibration/cameraParams.mat'
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
edgearr = [leftedgearr; rightedgearr];
figure()
imshow(edge)
%% Angle calculation
%Left edge angle
left1 = edgearr(edgearr(:,1) < 1460 &...
                edgearr(:,1) > 1200 &...
                edgearr(:,2) < 920 &...
                edgearr(:,2) > 740, :);
left2 = edgearr(edgearr(:,1) < 1740 &...
                edgearr(:,1) > 1730 &...
                edgearr(:,2) < 1800 &...
                edgearr(:,2) > 1720, :);

right1 = edgearr(edgearr(:,1) < 1520 &...
                edgearr(:,1) > 1280 &...
                edgearr(:,2) < 750 &...
                edgearr(:,2) > 580, :);
right2 = edgearr(edgearr(:,1) < 1918 &...
                edgearr(:,1) > 1913 &...
                edgearr(:,2) < 1800 &...
                edgearr(:,2) > 1720, :);

figure;
plot(edgearr(:,2),edgearr(:,1),'o', 'MarkerSize', 1)
hold on;
plot(left1(:,2), left1(:,1),'LineWidth',1.5)
plot(left2(:,2), left2(:,1),'LineWidth',1.5)
plot(right1(:,2), right1(:,1),'LineWidth',1.5)
plot(right2(:,2), right2(:,1),'LineWidth',1.5)
legend('p','L1','L2', 'R1', 'R2')
axis equal
                   
pL1 = polyfit(left1(:,2), left1(:,1), 1);
pL2 = polyfit(left2(:,2), left2(:,1), 1);
angleL1 = atand(pL1(1));
angleL2 = atand(pL2(1))
angleL = 180-angleL1+angleL2

pR1 = polyfit(right1(:,2), right1(:,1), 1);
pR2 = polyfit(right2(:,2), right2(:,1), 1);
angleR1 = atand(pR1(1));
angleR2 = atand(pR2(1))
angleR = 180-angleR1+angleR2


