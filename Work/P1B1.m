% Chris Woodruff 6/25/19
% This program runs the complete algorithm for calculating the bend angle
% by binarizing the image, identifying the centerline, and fitting a line
% to the centerline, taking the atan of the slope as the bend angle. 
% P1B1 107.26
% P1B2 125.23
%% Load Images
clc, clear all, close all;
% load image of tube bend
im = imread('white.png');
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
figure()
imshow(edge)
%% Angle calculation
%Left edge angle
left1 = leftedgearr(leftedgearr(:,1) < 1350 &...
                   leftedgearr(:,1) > 1110 &...
                   leftedgearr(:,2) < 680 &...
                   leftedgearr(:,2) > 620, :);
left2 = leftedgearr(leftedgearr(:,1) < 1758 &...
                   leftedgearr(:,1) > 1748 &...
                   leftedgearr(:,2) < 1877 &...
                   leftedgearr(:,2) > 1654, :);

figure;
plot(leftedgearr(:,2),leftedgearr(:,1),'o', 'MarkerSize', 1)
hold on;
plot(left1(:,2), left1(:,1),'LineWidth',1.5)
plot(left2(:,2), left2(:,1),'LineWidth',1.5)
legend('p','L1','L2')
axis equal
                   
pL1 = polyfit(left1(:,2), left1(:,1), 1);
pL2 = polyfit(left2(:,2), left2(:,1), 1);
angleL1 = atand(pL1(1))
angleL2 = atand(pL2(1))
180-angleL1+angleL2

