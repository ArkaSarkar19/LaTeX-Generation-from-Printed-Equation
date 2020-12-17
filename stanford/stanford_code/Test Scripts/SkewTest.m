% This script tests the fn_deskew2 function by reading in the "clean"
% equation images in the "LaTeX Equations" directory. Each image is named
% "eq[i]_hr.jpg" where i  is an integer. Currently, there is no eq4_hr.jpg.
% The script manually rotates each image and then runs fn_deskew2 to detect
% the rotation. fn_deskew2 outputs the deskewing angle, which is the
% negative of the angle used to skew the image. This script plots the
% absolute angular error in degrees.
addpath('../');
%% Read in the images
imgs = cell(13);
j = 1;
% equations 12 through 14 are similar to eq11. All 4 equations achieve poor
% results due to brevity of equation and prevalence of non horizontal
% lines.
% for i = 1:14
for i = 1:10
    if i == 4, continue, end;
    imgs{j} = rgb2gray(imread(['../Equations/Clean/eq' num2str(i) '_hr.jpg']));
    j = j+1;
end

%% Perform the tests
detected = zeros(size(-25:0.5:25)); % Holds the detected deskewing angle.
err = zeros(13,length(detected)); % Holds the absolute errors
for j = 1:13
    img = imgs{j};
    i = 1;
    for theta = -25:0.5:25
        
        % Rotate each image with bilinear interpretation and with
        % background filling.
        rotated = imrotate(img,theta);
        Mrot = ~imrotate(true(size(img)),theta,'bilinear');
        rotated(Mrot&~imclearborder(Mrot)) = 255;
        
        % Binarize using Otsu's method
        thresh = graythresh(rotated);
        bw_skew = im2bw(rotated, thresh);
        
        % Run algorithm
        [~,detected(i)] = fn_deskew2(bw_skew);
        if theta ==0
            err(j,i) = abs((detected(i)+theta));
        else
            err(j,i) = abs((detected(i)+theta));
        end
        i = i+1;
    end
end;
plot(-25:0.5:25,err(1:10,:))
legend('eq1','eq2','eq3','eq5','eq6','eq7','eq8','eq9','eq10')
title('Absolute Rotation Error')