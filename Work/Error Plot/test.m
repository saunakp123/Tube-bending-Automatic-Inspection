% P1B1 107.26
% P1B2 125.23 108.12
% P2B1 90.15
% P2B2 120.37 90.21
% P2B3 129.92 120.23 90.01
% P2B4 129.53
clc, clear all, close all;

name = 'P2B4'
% 1st value is the true angle
% 2nd - 5th value are the separation point for inner edge
% 6th - 9th value are the separation point for outer edge
% 10th value is the offset of the line to separate two edges
P1B1 = [107.26, 1092 603 1451 720, 1144 434 1541 568, 200];
P1B2 = [125.23, 1034 648 1595 1033, 1095 481 1732 921, -200];
P2B1 = [90.15, 229 630 1122 640, 228 455 1221 470, 500];
P2B2 = [120.37, 772 304 1559 788, 908 183 1637 624, -200];
P2B3 = [129.92, 668 1 1571 753, 949 1 1704 634, -600];
P2B4 = [129.53, 917 369 1623 1001, 1101 312 1744 869, -500];

tube = eval(name);
true_angle = tube(1);
leftind = tube(2:5);
rightind = tube(6:end-1);

%% Load Images
% load image of tube bend
im = imread(['./Images/Previous/', name, '.png']);
% load 'zero' image with no tube for background comparison
imz = imread('./Images/Previous/white.png');
% load saved camera parameters
% load './Previous work/Camera Calibration/cameraParams.mat'
% figure; imshow(im);


%% Image Filtering
% [white_u,newOrigin] = undistortImage(im,cameraParams);
% [white_uz,newOrigin] = undistortImage(imz,cameraParams);
diff = imz - im;
im2 = imbinarize(diff);
% figure; imshow(im2);
bw2=imfill(im2,'holes');
% figure;imshow(bw2);

edges = edge(bw2);
figure; imshow(edges);


%% Edge detection
% Define the separating line
k = tand(true_angle - 90);
b = tube(end);
x = 0:0.1:2000;
y = k*x + b;

% Separating the left & right edge indices
[i, j, v] = find(edges);
ind = [i j];
left = ind(ind(:, 1) >= leftind(2)  &...
            ind(:, 1) <= leftind(4) &...
            ind(:, 2) >= leftind(1) &...
            ind(:, 2) <= leftind(3) &...
            k*ind(:, 2) + b < ind(:,1), :);

right = ind(ind(:, 1) >= rightind(2) &...
            ind(:, 1) <= rightind(4) &...
            ind(:, 2) >= rightind(1) &...
            ind(:, 2) <= rightind(3) &...
            k*ind(:, 2) + b > ind(:,1), :);


figure;
plot(ind(:,2), ind(:,1),'o', 'MarkerSize', 1)
hold on;

plot(left(:,2), left(:,1),'LineWidth',1.5)
plot(right(:,2), right(:,1),'LineWidth',1.5)
plot(x, y, 'LineWidth',1.5)
legend('p','L','R', 'separate line')
axis equal


%% Angle calculation
% Start from the middle or edge
[lens, angle_errors, fit_errors] = getErrorFromMiddle(left, true_angle);
y = 0.1 * ones([1, length(lens)]);
figure('position', [0, 0, 1920, 1080])
subplot(2, 1, 1)
plot(lens, angle_errors);
hold on;
plot(lens, fit_errors);
plot(lens, y)
ylim([0, 0.5])
xlabel('Length', 'FontSize', 14)
ylabel('Error', 'FontSize', 14)
legend({'angle error', 'fit error', 'error limit'}, 'FontSize', 14)
title('Inner Edge', 'FontSize', 18)

[lens, angle_errors, fit_errors] = getError(right, true_angle);
y = 0.1 * ones([1, length(lens)]);
subplot(2, 1, 2)
plot(lens, angle_errors);
hold on;
plot(lens, fit_errors);
plot(lens, y);
ylim([0, 0.5])
xlabel('Length', 'FontSize', 14)
ylabel('Error', 'FontSize', 14)
legend({'angle error', 'fit error', 'error limit'}, 'FontSize', 14)
title('Outer Edge', 'FontSize', 18)


saveas(gcf, ['./Error Plot/35inch/', name, '.png'])










