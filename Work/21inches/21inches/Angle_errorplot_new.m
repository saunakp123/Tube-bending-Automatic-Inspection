% This program runs the complete algorithm for calculating the bend angle
% by binarizing the image, identifying the centerline, and fitting a line
% to the centerline, taking the atan of the slope as the bend angle. 
%% Load Images
clc, clear all, close all;
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
% undistort image and 'zero' image using calibrated camera parameters
[white_u,newOrigin] = undistortImage(im,cameraParams);
[white_uz,newOrigin] = undistortImage(imz,cameraParams);
% subtract 'zero' image from original
im=white_uz-white_u;
%im = imz - im ;
% figure;
% imshow(im);
% binarize image
im2=imbinarize(im);
im3 = imdilate(im2,strel('disk',5));
%im4 = bwmorph(im3,'thin',inf);
%figure; imshow(im2);
%figure; imshow(im3);
%figure; imshow(im4); 
% fill in any holes caused by noisy image
bw2=imfill(im3,'holes');
%figure;imshow(bw2);
%% Edge detection (stalled)
% % identifying edges 
% k=0;
% l=0;
% leftedge = 0*bw2 ;
% rightedge = leftedge ;
% for i=764:n_x 
%     for j = 848:n_y-1
%         if bw2(i,j)==0 && bw2(i,j+1)==1
%             k = k+1 ;
%             leftedge(i,j+1) = bw2(i,j+1) ;
%             leftedgearr(k,:) = [j+1,i] ;
%         end
%         if bw2(i,j)==1 && bw2(i,j+1)==0
%             l = l+1 ;
%             rightedge(i,j) = bw2(i,j) ;
%             rightedgearr(l,:) = [j,i] ;
%             continue ;
%         end
%     end
% end
% %imshow(rightedge) 
% %figure() ;
% %imshow(leftedge)
% edge = leftedge + rightedge ;
% figure()
% imshow(edge)
% edgearr = [leftedgearr; rightedgearr] ;
% %plot(edgearr(:,2), edgearr(:,1),'o')
%% Edge detection (stalled)
% identifying edges 
k=0;
l=0;
leftedge = 0*bw2 ;
rightedge = leftedge ;
[n_x n_y] = size(bw2);
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
edge = leftedge + rightedge ;
figure()
imshow(edge)
%% Edge detection
% low = 0.2 ;
% high = 0.6 ;
% sigma = sqrt(10) ;
% edge = edge(bw2,'canny',[low high],sigma);
% figure ;imshow(edge);
%% Angle calculation
%Storing array indices as left and right (X,Y) points 
Yl = leftedgearr(:,2);
Xl = leftedgearr(:,1);
Yr = rightedgearr(:,2);
Xr = rightedgearr(:,1);
%Plotting the Y's and X's
figure ;
plot(Yl,Xl);
axis equal
hold on
plot(Yr,Xr) ;
%
%For P1B1
%left_slant_llim = find(Yl==383) ; 

%For P1B2,P2b1,P2B3,P2B4(21) & P2B2(33.5)
left_slant_llim = find(Yl==50) ;

%For P2B3(33.5)
%left_slant_llim = find(Yl==280) ;

%For P1B2(33.5)
% left_slant_llim = find(Yl==855) ;
% right_slant_llim = find(Yr==750) ;
% %For P1B2(hor-33.5)
%  left_hor_llim = find(Yl==1710);
%  left_hor_ulim = find(Yl==1760);
% % 
%  right_hor_llim = find(Yr==1800);
%  right_hor_ulim = find(Yr==1850);

%For P1
 left_hor_llim = find(Yl==1830);
 left_hor_ulim = find(Yl==1880);
% 
 right_hor_llim = find(Yr==1560);
 right_hor_ulim = find(Yr==1610);

% %For P2B2
 left_hor_llim = find(Yl==1620);
 left_hor_ulim = find(Yl==1670);
% 
 right_hor_llim = find(Yr==1785);
 right_hor_ulim = find(Yr==1835);
% 
%For P2B3
left_hor_llim = find(Yl==1525);
left_hor_ulim = find(Yl==1575);
% %For P2B3(33.5)
% left_hor_llim = find(Yl==1780);
% left_hor_ulim = find(Yl==1830);
% 
right_hor_llim = find(Yr==1860);
right_hor_ulim = find(Yr==1910);

For P2B4
left_hor_llim = find(Yl==1880);
left_hor_ulim = find(Yl==1930);

right_hor_llim = find(Yr==1930);
right_hor_ulim = find(Yr==1980);

%For P2B4(33.5)

% left_slant_llim = find(Yl==660) ;
% right_slant_llim = find(Yr==530) ;
% 
% left_hor_llim = find(Yl==1660);
% left_hor_ulim = find(Yl==1710);
% 
% right_hor_llim = find(Yr==1650);
% right_hor_ulim = find(Yr==1700);
j = 1 ;
start_p = 0 ;
%start_p = [0:25:700] ;
%for j = 1:length(start_p)
gp_left_slant = [start_p(j)+left_slant_llim start_p(j)+left_slant_llim+20] ; %good points left slant
d_left = abs(gp_left_slant(1,1) - gp_left_slant(1,2)) ;
gp_left_hor = [left_hor_llim left_hor_ulim] ; %good points left horizontal
gp_right_slant = [start_p(j)+Yr(10) start_p(j)+Yr(20)] ; %good points right slant
d_right = abs(gp_right_slant(1,1) - gp_right_slant(1,2)) ;
gp_right_hor = [right_hor_llim right_hor_ulim] ; %good points right horizontal
d = 0:10:1000 ;
n = length(d) ;
% Adding points and increasing the length of the strip for angle calculation
for i = 1:n
    Yl_1 = Yl(gp_left_slant(1,1):gp_left_slant(1,2)+0.5*d(i));
    Xl_1 = Xl(gp_left_slant(1,1):gp_left_slant(1,2)+0.5*d(i));
    Yl_2 = Yl(gp_left_hor(1,1):gp_left_hor(1,2));
    Xl_2 = Xl(gp_left_hor(1,1):gp_left_hor(1,2));
    Yr_1 = Yr(gp_right_slant(1,1):gp_right_slant(1,2)+0.5*d(i));
    Xr_1 = Xr(gp_right_slant(1,1):gp_right_slant(1,2)+0.5*d(i));
    Yr_2 = Yr(gp_right_hor(1,1):gp_right_hor(1,2));
    Xr_2 = Xr(gp_right_hor(1,1):gp_right_hor(1,2));
    L_left(i) = d_left + 0.5*d(i) ;
    L_right(i) = d_right + 0.5*d(i) ;
    %Curvefitting to get slope of left edge line(before smoothing)
    p = polyfit(Yl_1,Xl_1,1) ;
    p1 = polyfit(Yl_2,Xl_2,1) ; % using right edge horizontal points 
    angle_l = atan(p(1))*180/pi ;
    angle1_l = atan(p1(1))*180/pi ;
    left_angle(i) = 180-angle_l+angle1_l ;
    error_left_edge(i) = abs(left_angle(i)-actual) ;    
    %Curvefitting to get slope of right edge line(before smoothing)
    p = polyfit(Yr_1,Xr_1,1) ;
    p1 = polyfit(Yr_2,Xr_2,1) ;
    angle_r = atan(p(1))*180/pi ;
    angle1_r = atan(p1(1))*180/pi ;
    right_angle(i) = 180-angle_r+angle1_r ;
    error_right_edge(i) = abs(right_angle(i)-actual) ;     
end
error_right_edge_index = error_right_edge<0.5 ;
error_left_edge_index = error_left_edge<0.5 ;
error_right_edge=error_right_edge(error_right_edge_index)
error_left_edge = error_left_edge(error_left_edge_index)
L_left = L_left(error_left_edge_index) ;
L_right = L_right(error_right_edge_index);
%% Error plotting
[elmin,li] = min(error_left_edge) ;
[ermin,ri] = min(error_right_edge) ;
figure ;
plot(L_left(li),error_left_edge(li),'ro','Markersize',10) ;
hold on ;
plot(L_right(ri),error_right_edge(ri),'bo','Markersize',10) ;
plot(L_left,error_left_edge,'b-') ; 
plot(L_right,error_right_edge,'r-') ;
xlabel('Length') ;
ylabel('Error') ;
title(sprintf('Error plot: %s(%0.3f in), inner start: %d & outer start: %d',upper(part),camera_dist,start_p(j)+left_slant_llim,start_p(j)+Yr(10))) ;
legend(sprintf('Inner strip minimum error = %0.3f & %d pixels',elmin,L_left(li)),sprintf('Outer strip minimum error  = %0.3f & %d pixels',ermin,L_right(ri)),'Inner strip','Outer strip');
legend('Location', 'southoutside') ;
savefig(sprintf('%s_figure %d',upper(part),j)) ;
hold off
%end