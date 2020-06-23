% P1B1 107.26
% P1B2 125.23 108.12
% P2B1 90.15
% P2B2 120.37 90.21
% P2B3 129.92 120.23 90.01
% P2B4 129.53
clc, clear, close all;

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
warning('off')
%Choose a number to get an image 
part = input('Enter the part name: ','s') ;
camera_dist = input('Enter the camera distance: ') ;
[im,imz,actual,load_path] = chooseimage(part,camera_dist) ;
load(load_path) ;
% display raw image
 figure;
 imshow(im);
%% Image Filtering
 [im,newOrigin] = undistortImage(im,cameraParams);
 [imz,newOrigin] = undistortImage(imz,cameraParams);
diff = imz - im;
im2 = imbinarize(diff);
im3 = imdilate(im2,strel('disk',5)) ;
% figure; imshow(im2);
bw2=imfill(im3,'holes');
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

%Inner 
%figure() ;
for i = 1:10:(length(right)/2)-2
    p1 = right(i,:) ; % Point 1 (From the beginning)
    for j = length(right):-10:(length(right)/2)+2    
        p2 = right(j,:) ; % Point 2 (From the end)
        p = polyfit([p1(2) p2(2)],[p1(1) p2(1)],1) ; % Finding the slope of the 2 points
        dis(i,length(right)+1-j) = norm(right(i,:)-right(j,:)) ;
        ang = atand(p(1)) ; % Finding the angle of the line from the 2 points 
        ang = ang+90 ; % Adding 90 degrees
        err(i,length(right)+1-j) = (ang-true_angle) ; % Angle error for each starting point and ending points
        %plot([p1(2) p2(2)],[p1(1) p2(1)]) ;
        %hold on ;
    end
end
dis2 = dis(:) ;
[dis3,ind] = sort(dis2) ;
err2 = err(:) ;
err3 = err2(ind) ;
y = 0.1*ones(size(dis2));
figure ;
plot(dis3,abs(err3)) ;
hold on
plot(dis2,y,'r--')
for i = 1:10:(length(left)/2)-2
    p1 = left(i,:) ; % Point 1 (From the beginning)
    for j = length(left):-10:(length(left)/2)+2    
        p2 = left(j,:) ; % Point 2 (From the end)
        p = polyfit([p1(2) p2(2)],[p1(1) p2(1)],1) ; % Finding the slope of the 2 points
        dis(i,length(left)+1-j) = norm(left(i,:)-left(j,:)) ;
        ang = atand(p(1)) ; % Finding the angle of the line from the 2 points 
        ang = ang+90 ; % Adding 90 degrees
        err(i,length(left)+1-j) = (ang-true_angle) ; % Angle error for each starting point and ending points
        %plot([p1(2) p2(2)],[p1(1) p2(1)]) ;
        %hold on ;
    end
end
dis2 = dis(:) ;
[dis3,ind] = sort(dis2) ;
err2 = err(:) ;
err3 = err2(ind) ;
 plot(dis3,abs(err3)) ;
title(sprintf('Angle error vs Pixels between points: %s',name))
legend('Error from line fit(outer)','Threshold error','Error from line fit(inner)')
ylabel('Angle error')
xlabel('Pixels between points')