function [ identifier ] = create_indent(image)

showFigs = true;
% Parameters
k = 8;
% IN2, R1-R7, D2-D8, 6 Hu Moments
identifier = zeros(1,2*k+6);

inverse_image = ones(size(image)) - image;

figure, imshow(image, []);
title("image")

figure, imshow(inverse_image, []);
title("inverse_image")
% If white border on any edge, pad with more 1s
pad = 3;
% top
if(sum(inverse_image(1,:)) == 0)
    image = [ones(pad,size(image,2)) ;image];
end
% bottom
if(sum(inverse_image(end,:)) == 0)
    image = [image; ones(pad,size(image,2))];
end
% left
if(sum(inverse_image(:,1)) == 0)
    image = [ones(size(image,1),pad) image];
end
% right
if(sum(inverse_image(:,end)) == 0)
    image = [image ones(size(image,1),pad)];
end
% Recalculate for centroid
inverse_image = ones(size(image)) - image;

% Total number of pixels in character
N = sum(~image(:));
% Find x,y positions of 0 elements (character)
[y, x] = find(~image);

% Find centroid
centroid = regionprops(inverse_image,'centroid');
centroid = cat(1, centroid.Centroid);
centroid_x = centroid(1);
centroid_y = centroid(2);

moment_of_inertia = sum((x - centroid_x).^2 + (y - centroid_y).^2) /sum(~image(:))^2;;
identifier(1) = moment_of_inertia;

hu_moments =  calculate_hu_moments(inverse_image);
identifier(2*k+1:end) = hu_moments(2:end);


y = size(image,1);
x = size(image,2);
dr = max(max(x - centroid_x,centroid_x), max(y - centroid_y,centroid_y)) / (k+1);
[c, r] = meshgrid(1:x, 1:y);


coding(k).circVec=[];
for i = 1:k
    rad = dr * i;
    % Create circle line of logical 1s to extract template
    rad1 = sqrt((c-centroid_x).^2+(r-centroid_y).^2);
    rad2 = sqrt((c-centroid_x).^2+(r-centroid_y).^2);
    Circle = xor(rad1<=rad, rad2<=(rad-1));
    % Extract circular template linear indices (sort based on theta)
    cidx = find(Circle);
    % Create matrix linear indices, corresponding x and y values
    x = c(cidx)-centroid_x;
    y = r(cidx)-centroid_y;
    vals = [cidx x y zeros(size(cidx,1),1)];
    % Add corresponding theta values
    [thetha, rho] = cart2pol(vals(:,2), vals(:,3));
    vals(:,4) = thetha;
    % Sort based on theta values so in order around circle
    [~, order] = sort(vals(:,4));
    sortedvals = vals(order,:);
    % Use linear indices to extract desired values in correct theta order
    circVec = image(sortedvals(:,1));
    
    % Show circles overlaid on character
%     if(showFigs)
%         tempcirc(Circle) = .5;
%     end
    
    % Extract values from character based on circular indices. Start at
    %   0deg going CCW.
    circVec = imopen(circVec,ones(2,1));
    
    % Save vector for debugging or future use
    coding(i).circ = circVec;
    
    % If character is so small nothing is extracted, leave ident as 0
    if(~isempty(circVec))
        % Count number of sections of at least 2 0s. This is the number of
        % character sections the circle goes through.
        arr = [1 1 circVec'];
        cnt = getcount(arr);
        coding(i).count = length(cnt(diff([1 cnt])~=1));

        % Save count to identifier vector
        identifier(1+i) = coding(i).count;
        if(showFigs)
            subplot(2,4,i);
            imshow(coding(i).circ);
            str = sprintf('%d',coding(i).count);
            title(str);
        end

        % Find 2 longest arcs of background and divide diff by total length
        % If 1s (background) wrap around, move all 1s to one side.
        circ = length(circVec);
        if(circVec(1) == 1 && circVec(end) == 1)
            idx = find(circVec==0,1,'last');
            circVec = [circVec(idx+1:end); circVec(1:idx)];
        end
        % Vector of length of each section of 1s
        B = [0 circVec' 0]; % Pad with 0s for diff
        bgrd_len = find(diff(B)==-1) - find(diff(B)==1);
        % If no crossings, set ratio to 0
        if(coding(i).count < 1);
            d2 = 0;
            d1 = 0;
        else
            % Get two longest sections (arcs)
            [d2, idx_d2] = max(bgrd_len);
            bgrd_len(idx_d2)=NaN;
            d1 = max(bgrd_len);
        end
        % Find ratio of difference of two longest arcs by circumference
        coding(i).ratio = (d2-d1) / circ;
        if(i > 1) % Only keep k-1 largest ratios
            if(~isnan(coding(i).ratio))
                identifier(k+i) = coding(i).ratio;
            end
        end

    end
end


end




function hu_arr = calculate_hu_moments(inverse_image)
eta_mat = SI_Moment(inverse_image) ;
hu_arr = Hu_Moments(eta_mat);
end

function count_zero = getcount(arr)

l = length(arr);
ctr = 1;
count_zero = zeros(l);
for i=1:l-1
    if(arr(i) == 0 && arr(i+1) == 0)
        count_zero(ctr) = i;
        ctr = ctr + 1;
    end
end
count_zero = count_zero(1:ctr-1);
end



%% 
% Author:   Vishnu Muralidharan
% University of Alabama in Huntsville
% Copyright (c) 2015, Vishnu Muralidharan
% All rights reserved.

function eta = SI_Moment(image, mask)
% Function to calculate the scale-invariant moment of interested image region
% Author:   Vishnu Muralidharan
% University of Alabama in Huntsville
% Copyright (c) 2015, Vishnu Muralidharan
% All rights reserved.


% Inputs:   image: image: input image for which moments need to be calculated
%           mask: specifying this allows you to calculate moments for a
%           specified region
%           
% Outputs:  cen_mmt = central moment of the specifed order fot the image
% Reference:  Visual Pattern Recognition by Moment Invariants


image = double(image);
if ~exist('mask','var')
    mask = ones(size(image,1),size(image,2)); %if mask is not defined select the whole image
end

% computation of central moments upto order order 3
for i=1:1:4
    for j=1:1:4
        nu(i,j) = Centr_Moment(image, mask,i-1,j-1);
    end
end

% computation of scale invariant moments using central monets of upto order
% 3
eta = zeros(3,3);
for i=1:1:4
    for j=1:1:4
        if i+j >= 4
            eta(i,j) = (double(nu(i,j))/(double(nu(1,1)).^(double((i+j)/2)))); %scale invariant moment matrix
        end
    end
end
end


function inv_moments = Hu_Moments(eta)
% Function to calculate the Hu's moments of interested image region
% Author:   Vishnu Muralidharan
% University of Alabama in Huntsville
% Copyright (c) 2015, Vishnu Muralidharan
% All rights reserved.


% Inputs:   eta: scale-invariant moment matrix of order 3
% Outputs:  inv_moments = array containing 7 invariant Hu's moments
% Reference:  Visual Pattern Recognition by Moment Invariants


%Calculation of various invariant Hu's moments
inv_moments(1) = eta(3,1) + eta(1,3);
inv_moments(2) = (eta(3,1) - eta(1,3))^2 + (4*eta(2,2)^2);
inv_moments(3) = (eta(4,1) - 3*eta(2,3))^2 + (3*eta(3,2) - eta(1,4))^2;
inv_moments(4) = (eta(4,1) + eta(2,3))^2 + (eta(3,1) + eta(1,4))^2;
inv_moments(5) = (eta(4,1) - 3*eta(2,3))*(eta(4,1) + eta(2,3))*((eta(4,1) + eta(2,3))^2 - 3*((eta(3,2) + eta(1,4))^2)) + (3*(eta(3,2) - eta(1,4)))*(eta(3,2) + eta(1,4))*(3*(eta(4,1) + eta(2,3))^2 - (eta(3,2) + eta(1,4))^2);
inv_moments(6) = (eta(3,1) - eta(1,3))*((eta(4,1)+eta(2,3))^2 - (eta(3,2)+ eta(1,4))^2) + 4*eta(2,2)*((eta(4,1) + eta(2,3))*(eta(3,2) + eta(1,4)));
inv_moments(7) = (3*eta(3,2) - eta(1,4))*(eta(4,1) + eta(2,3))*((eta(4,1) + eta(2,3))^2 - 3*(eta(3,2)-eta(1,4))^2) - (eta(4,1) - 3*eta(2,3))*(eta(3,2) + eta(1,4))*(3*(eta(4,1) + eta(2,3))^2 - (eta(3,2) + eta(1,4))^2);

end


function cen_mmt = Centr_Moment(image,mask,p,q)
% Function to calculate the central moment of interested image region
% Author:   Vishnu Muralidharan
% University of Alabama in Huntsville
% Copyright (c) 2015, Vishnu Muralidharan
% All rights reserved.


% Inputs:   image: image: input image for which moments need to be calculated
%           mask: specifying this allows you to calculate moments for a
%           specified region
%           p,q: order of moments to be calculated
% Outputs:  cen_mmt = central moment of the specifed order fot the image
% Reference:  Visual Pattern Recognition by Moment Invariants


if ~exist('mask','var')
    mask = ones(size(image,1),size(image,2)); %if mask is not spcified, select the whole image
end

image = double(image);

%moments necessary to compute components of centroid
m10 = moment(image,mask,1,0); 
m01 = moment(image,mask,0,1);
m00 = moment(image,mask,0,0);

%components of centroid
x_cen = floor(m10/m00);
y_cen = floor(m01/m00);

cen_mmt =0;

for i=1:1:size(mask,1)
    for j=1:1:size(mask,2)
        if mask(i,j) == 1
            cen_mmt = cen_mmt + (double(image(i,j))*((i-x_cen)^p)*((j-y_cen)^q)); %calculating central moment
        end
    end
end
end


function m = moment(image,mask,p,q)
% Function to calculate any ordinary moment of the intersted image region
% Author:   Vishnu Muralidharan
% University of Alabama in Huntsville
% Copyright (c) 2015, Vishnu Muralidharan
% All rights reserved.


% Inputs:   image: input image for which moments need to be calculated
%           mask: specifying this allows you to calculate moments for a
%           specified region
%           p,q: order of moments to be calculated
% Outputs:  m = moment of the specifed order fot the image
% Reference:  Visual Pattern Recognition by Moment Invariants


if ~exist('mask','var')
    mask = ones(size(image,1),size(image,2));   %if mask is not specified, select the whole image
end

image = double(image);
m=0; 
for i=1:1:size(mask,1)
    for j=1:1:size(mask,2)
        if mask(i,j) == 1
            m = m + (double((image(i,j))*(i^p)*(j^q))) ; %moment calculation
        end
    end
end
end