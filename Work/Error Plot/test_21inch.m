% P1B1 107.26
% P1B2 125.23 108.12
% P2B1 90.15
% P2B2 120.37 90.21
% P2B3 129.92 120.23 90.01
% P2B4 129.53
clc, clear all, close all;

name = 'P1B1';
% 1st value is the true angle
% 2nd - 5th value are the separation point for inner edge
% 6th - 9th value are the separation point for outer edge
% 10th value is the offset of the line to separate two edges
P1B1 = [107.26, 561 892 1104 1058, 659 607 1125 752, 500];
P1B2 = [125.23, 665 845 1351 1330, 728 537 1585 1143, 200];
P2B2 = [120.37, 190 570 1259 1193, 385 336 1369 911, 200];
P2B3 = [129.92, 1 356 1338 1461, 43 1 1545 1243, 200];
P2B4 = [129.53, 508 812 1287 1448, 611 502 1492 1225, 200];

tube = eval(name);
true_angle = tube(1);
leftind = tube(2:5);
rightind = tube(6:end-1);

%% Load Images
% load image of tube bend
im = imread(['./Images/21inch/', name, '.png']);
% load 'zero' image with no tube for background comparison
imz = imread('./Images/21inch/white.png');
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

[i, j, v] = find(edges);
ind = [i j];
% Separating the left & right edge indices
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
% [lens, angle_errors, fit_errors] = getErrorFromEdge(left, true_angle);
y = 0.1 * ones([1, length(lens)]);
max_length = length(lens);
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

[lens, angle_errors, fit_errors] = getErrorFromEdge(right, true_angle);
y = 0.1 * ones([1, length(lens)]);
max_length = max(max_length, length(lens));
subplot(2, 1, 2)
plot(lens, angle_errors);
hold on;
plot(lens, fit_errors);
plot(lens, y);
ylim([0, 0.5])
xlabel('Length', 'FontSize', 14)
ylabel('Error', 'FontSize', 14)
% legend({'angle error', 'fit error', 'error limit'}, 'FontSize', 14)
title('Outer Edge', 'FontSize', 18)


% saveas(gcf, ['./Error Plot/21inch/', name, '.png'])


%% Idealized Error
angle = true_angle - 90;

l_choice = 10:20:max_length;

x = 1:2000;
y = tand(angle) * x;
delta_y = 0:0.01:1;
y = y + reshape(delta_y, [], 1);
y = floor(y);

fit_angles = zeros(length(delta_y), length(l_choice));
for ln = 1:length(l_choice)
    l = l_choice(ln);
    xcut = round(l * cosd(angle));
    xp = x(1:xcut);
    yps = y(:, 1:xcut);
    for i = 1:size(yps, 1)
        yp = yps(i, :);
        za = polyfit(xp, yp, 1);
        fit_angles(i, ln) = atand(za(1));
    end
end

errors = abs(fit_angles - angle);
error_mean = mean(errors, 1);
error_max = max(errors, [], 1);
error_min = min(errors, [], 1);

plot(l_choice, error_max)
plot(l_choice, error_min)



