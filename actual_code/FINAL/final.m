%% 
close all;
clear
clc
fileName = 'wrong.jpg';
Equation_image = imread(fileName);
figure("Name","Original Image");                                                              
imshow(Equation_image);
%% 
load('red_charPalette_withText_demo2.mat');
load('red_charPalette_Classifier_demo2.mat');
%% 
close all; 


bin_image=binarization(Equation_image);
figure("Name","Binary Image");                                                              
imshow(bin_image);
%% 
close all; 
[Corrected_image, ~] = skew_correction(bin_image);



%% 

close all;
[ characters_image,characters_centroids,characters_boxes ] = segmentation(Corrected_image);
characters_ident(size(characters_image,1)).ident=[];
%% 
close all;

for i = 1:length(characters_image)
   character_ident(i).ident = Find_characteristics(characters_image(i).img); 
end
%% 

close all;
for i = 1:length(characters_image)

    
    idx_matched = similarity_function(character_ident(i).ident, "Manhattan");
    character_ident(i).char = chars(X_orig(idx_matched,end)).char;
     
 
  
    
    if(true)
        figure("Name",strcat("character",int2str(i)))
        imshow(characters_image(i).img);
        title('Input');
        
       figure("Name",strcat("Matched character",int2str(i)))
        title('Matched');

        imshow(chars(X_orig(idx_matched,end)).img);
                title('Matched');

       
    end
end




%% 
close all;
blank.val=character_ident;
clc
final_string = Form_equation(blank,characters_boxes,characters_centroids);
%% 
final_string
figure("Name","Original Image");                                                              
imshow(Equation_image);

%% 



function binary_image = binarization(img)
gray_image=rgb2gray(img);
thresh = graythresh(gray_image); %gray thresh returns the otsu threshold

[m,n]=size(gray_image);
binary_image=zeros(size(gray_image));

for i=1:m
    for j=1:n
        if (gray_image(i,j)>(255*thresh))
            binary_image(i,j)=1;
        end
    
    end

end
binary_image=logical(binary_image); %convert everything in foreground to 1 and background to 0
end
function [ corrected_img, angle_correction ] = skew_correction( bw_img )


BW = edge(bw_img,'canny');
figure('Name',"Canny edges");
imshow(BW);

[H,T,rho] = hough(BW);
P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));

lines = houghlines(BW,T,rho,P,'FillGap',5,'MinLength',7);

figure('Name','Hough lines'), imshow(bw_img), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red'); %%change color of the longest segment



orientations = P(:,2);
if length(unique(orientations)) ~= 4    
    final_orien = mode(P(:,2));
else
    final_orien = orientations(1);

end

    
angle_correction = T(final_orien)-90;

if (angle_correction<0)
    angle_correction=mod(angle_correction,-90);
else
    angle_correction=mod(angle_correction,90);
 
end


corrected_img = imrotate(bw_img,angle_correction,'bilinear');
Mrot = ~imrotate(true(size(bw_img)),angle_correction,'bilinear');%to fix the background.
corrected_img(Mrot&~imclearborder(Mrot)) = true;
if (abs(angle_correction) > 2)
figure("Name","Corrected Image")
imshow(corrected_img)
end


end
function [ characters_image,characters_centroids,characters_boxes ] = segmentation(img)





exp = imerode(ones(size(img)) - img, ones(3,3));

figure("Name","Eroded image"), imshow(exp, []);

img_edges = xor(imerode(ones(size(img)) - img, ones(3,3)), ones(size(img)) - img);



Ctr = regionprops(img_edges, 'Centroid');
BoundingB = regionprops(img_edges, 'BoundingBox');



centroids = cat(1, Ctr.Centroid);

boundingboxes = cat(1, BoundingB.BoundingBox);
boundingboxes = floor(boundingboxes);
boundingboxes(:,3:4) = boundingboxes(:,3:4) + 1;


characters_image(size(centroids,1)).img=[];

dummy=zeros(size(centroids,1),1);
for i =1:size(dummy,1)
    for j=1:size(dummy,1)
        if (i~=j)
            if ((boundingboxes(i,2)<=boundingboxes(j,2))&&(boundingboxes(i,2)+boundingboxes(i,4)>boundingboxes(j,2)+boundingboxes(j,4))&&(boundingboxes(i,1)<=boundingboxes(j,1))&&(boundingboxes(i,1)+boundingboxes(i,3)>boundingboxes(j,1)+boundingboxes(j,3)))
                dummy(j)=1;
            end

        end
    end
    

end
counter=0;
for i=1:size(dummy)
    if (dummy(i)==1)
        characters_image(i-counter)=[];
        boundingboxes(i-counter,:)=[];
        centroids(i-counter,:)=[];
        counter=counter+1;
    end
    
        
end



for i = 1 : size(centroids,1)
    
    characters_image(i).img = img(boundingboxes(i,2):boundingboxes(i,2)+boundingboxes(i,4),boundingboxes(i,1):boundingboxes(i,1)+boundingboxes(i,3),:);
    
    
end



figure("Name","Centroids");

imshow(img_edges);

for i = 1:size(centroids, 1) 
    hold on ;
    x = centroids(i,1);
    y = centroids(i,2);
    text(x, y, '*' ,'Color', 'red', 'FontSize', 10);
    hold off;

end
figure("Name","Bounding Boxes");

imshow(img_edges);

for i = 1:size(boundingboxes, 1) 
    hold on ;
    rectangle('position',boundingboxes(i,:),'Edgecolor','y')

    hold off;

end

characters_centroids=centroids;
characters_boxes=boundingboxes;
end
function [ identifier ] = Find_characteristics(image)

k = 8;
identifier = zeros(1,2*k+6);

inverse_image = ones(size(image)) - image;

figure, imshow(image, []);
title("image")

figure, imshow(inverse_image, []);
title("inverse_image")
pad = 3;
if(sum(inverse_image(1,:)) == 0)
    image = [ones(pad,size(image,2)) ;image];
end
if(sum(inverse_image(end,:)) == 0)
    image = [image; ones(pad,size(image,2))];
end
if(sum(inverse_image(:,1)) == 0)
    image = [ones(size(image,1),pad) image];
end
if(sum(inverse_image(:,end)) == 0)
    image = [image ones(size(image,1),pad)];
end
inverse_image = ones(size(image)) - image;
Topo_circle=image;
N = sum(~image(:));
[y, x] = find(~image);

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
    rad1 = sqrt((c-centroid_x).^2+(r-centroid_y).^2);
    rad2 = sqrt((c-centroid_x).^2+(r-centroid_y).^2);
    Circle = xor(rad1<=rad, rad2<=(rad-1));
    cidx = find(Circle);
    x = c(cidx)-centroid_x;
    y = r(cidx)-centroid_y;
    vals = [cidx x y zeros(size(cidx,1),1)];
    [thetha, rho] = cart2pol(vals(:,2), vals(:,3));
    vals(:,4) = thetha;
    [~, order] = sort(vals(:,4));
    sortedvals = vals(order,:);
    circVec = image(sortedvals(:,1));
    
   
     Topo_circle(Circle) = .5;
    
    circVec = imopen(circVec,ones(2,1));
    
    coding(i).circ = circVec;
    
    
    if(~isempty(circVec))
        
        arr = [1 1 circVec'];
        cnt = getcount(arr);
        coding(i).count = length(cnt(diff([1 cnt])~=1));

        
        identifier(1+i) = coding(i).count;
        

        
        circ = length(circVec);
        if(circVec(1) == 1 && circVec(end) == 1)
            idx = find(circVec==0,1,'last');
            circVec = [circVec(idx+1:end); circVec(1:idx)];
        end
        B = [0 circVec' 0]; 
        bgrd_len = find(diff(B)==-1) - find(diff(B)==1);
        if(coding(i).count < 1);
            d2 = 0;
            d1 = 0;
        else
            [d2, idx_d2] = max(bgrd_len);
            bgrd_len(idx_d2)=NaN;
            d1 = max(bgrd_len);
        end
        coding(i).ratio = (d2-d1) / circ;
        if(i > 1) 
            if(~isnan(coding(i).ratio))
                identifier(k+i) = coding(i).ratio;
            end
        end

    end
end
figure();
imshow(Topo_circle)


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

for i=1:1:4
    for j=1:1:4
        nu(i,j) = Centr_Moment(image, mask,i-1,j-1);
    end
end
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

function character = similarity_function(img_vector, type)

if(type == "Manhattan")
    
    Orig_vectors = load('red_charPalette_Classifier_demo2.mat').X_orig;
    character_labels = load('red_charPalette_withText_demo2.mat');
    
    [r,c] = size(Orig_vectors);
    
    min_distance = 9999;
    min_label = 0;
    
    for i=1:r
        arr = Orig_vectors(i,1:end-1);
        label = Orig_vectors(i,end);
        dist = sum(abs(arr - img_vector));
        if(dist < min_distance ) 
            min_distance = dist;
            min_label = label;
        end
    end
end


if(type == "Eucledian")
    
    Orig_vectors = load('red_charPalette_Classifier_demo2.mat').X_orig;
    character_labels = load('red_charPalette_withText_demo2.mat');
    
    [r,c] = size(Orig_vectors);
    
    min_distance = 9999;
    min_label = 0;
    
    for i=1:r
        arr = Orig_vectors(i,1:end-1);
        label = Orig_vectors(i,end);
        dist = sum((arr - img_vector).^2);
        if(dist < min_distance ) 
            min_distance = dist;
            min_label = label;
        end
    end
end

character = min_label;
end
function Final_eqn = Form_equation( characters,boundingboxes,centroids )

control = Get_symbols;
limit_control = {'int','sum','prod','lim'}; 

space = 10; 

Final_eqn = '';

chars = characters.val;
num_chars = length(chars);



boxes = zeros(4,num_chars);

for i = 1:num_chars
    boxes(:,i) = boundingboxes(i,:);
end
boxes;


[~, idxs] = sort(boxes(1,:));
chars = chars(idxs);
boxes = boxes(:,idxs);

chars = Preprocess(chars,boxes);

if length(chars) ~= num_chars
    num_chars = length(chars);
    boxes = zeros(4,num_chars);
    for i = 1:num_chars
        boxes(:,i) = chars(i).boundingbox;
    end
end

prev_centroid_y_coord = NaN;
prev_fraction = false; 
is_super = false; 
is_sub = false; 

i = 1;
while i <= num_chars

    detected = chars(i).char;
    
    top_left = boxes(1,i);
    top_right = top_left+boxes(3,i);
    

    
    overlap_indexs = boxes(1,:)>=top_left & boxes(1,:)<=top_right;
    cur_overlap_boxes = boxes(:,overlap_indexs);
    cur_overlap_top = min(cur_overlap_boxes(2,:));
    cur_overlap_bottom = max(cur_overlap_boxes(2,:)+cur_overlap_boxes(4,:));
    total_height = cur_overlap_bottom - cur_overlap_top;
    overlap_indexs = boxes(1,:)>=top_left & boxes(1,:)<=top_right &...
    boxes(4,:) < 0.7*total_height;
    

    cur_overlap_boxes = boxes(:,overlap_indexs);
    topleft_list = cur_overlap_boxes(1,:);
    topright_list = cur_overlap_boxes(1,:)+cur_overlap_boxes(3,:);
    top_list = cur_overlap_boxes(2,:);
    bottom_list = cur_overlap_boxes(2,:) + cur_overlap_boxes(4,:);
    
    
    while ~isempty(topleft_list) &&...
            (min(topleft_list) < top_left || max(topright_list) > top_right)
        top_left = min(topleft_list);
        top_right = max(topright_list);
        total_height = max(bottom_list)-min(top_list);
        overlap_indexs = boxes(1,:)>=top_left & boxes(1,:)<=top_right &...
            boxes(4,:) < 0.7*total_height;
        cur_overlap_boxes = boxes(:,overlap_indexs);
        topleft_list = cur_overlap_boxes(1,:);
        topright_list = cur_overlap_boxes(1,:)+cur_overlap_boxes(3,:);
    end
    
    

    overlap_indexs(i) = false;
    num_overlaps = sum(overlap_indexs);
    
    if (num_overlaps ==0)
        
            
            if i-1 > 0 && ~prev_fraction
                ll_corner = boxes(2,i)+boxes(4,i);
                ul_corner = boxes(2,i);
                
                if ll_corner <= ceil(prev_centroid_y_coord)
                    
                    if ~strcmp(detected,'-') || ll_corner < (prev_centroid_y_coord-5*boxes(4,i))
                        
                        is_super = true;
                        Final_eqn = [strtrim(Final_eqn) '^{'];
                        
                        
                        prev_centroid_y_coord = centroids(i,2);
                        
                        if i+1<= num_chars
                            dist_to_next = boxes(1,i+1)-top_right;
                        else
                            dist_to_next = NaN;
                        end
                        detected = help_get_char(detected,i,chars,dist_to_next);
                    end
                elseif ul_corner >= floor(prev_centroid_y_coord)
                    if is_super
                        is_super = false;
                        Final_eqn = [strtrim(Final_eqn) '}'];
                        prev_fraction = false;
                        
                        prev_centroid_y_coord = centroids(i,2);
                        
                        if i+1<= num_chars
                            dist_to_next = boxes(1,i+1)-top_right;
                        else
                            dist_to_next = NaN;
                        end
                        detected = help_get_char(detected,i,chars,dist_to_next);
                    else
                        is_sub = true;
                        Final_eqn = [strtrim(Final_eqn) '_{'];
                        
                        if i+1<= num_chars
                            dist_to_next = boxes(1,i+1)-top_right;
                        else
                            dist_to_next = NaN;
                        end
                        detected = help_get_char(detected,i,chars,dist_to_next);
                        
                        prev_centroid_y_coord = centroids(i,2);
                    end
                else
                    prev_fraction = false;
                    
                    prev_centroid_y_coord = centroids(i,2);
                    
                    if i+1<= num_chars
                        dist_to_next = boxes(1,i+1)-top_right;
                    else
                        dist_to_next = NaN;
                    end
                    detected = help_get_char(detected,i,chars,dist_to_next);
                end
            else
                prev_fraction = false;
                
                prev_centroid_y_coord = centroids(i,2);
                
                if i+1<= num_chars
                    dist_to_next = boxes(1,i+1)-top_right;
                else
                    dist_to_next = NaN;
                end
                detected = help_get_char(detected,i,chars,dist_to_next);
            end
            Final_eqn = [Final_eqn detected];
    elseif (num_overlaps==1)
            assert(overlap_indexs(i+1))
            overlap_ul = boxes(1,i+1);
            if abs(overlap_ul-top_right) <= 1
                if any(strcmp(detected,control))
                    detected = ['\' detected ' '];
                end
                
                if i+1 <= num_chars && help_is_letter(detected) && help_is_letter(chars(i+1).char)
                    dist_to_next = boxes(1,i+1)-top_right;
                    if dist_to_next >= space
                        detected = [detected '\,'];
                    end
                end
                Final_eqn = [Final_eqn detected];
                
                prev_centroid_y_coord = chars(i).centroid(2);
            else
                overlap_char = chars(i+1).char;
                
                if strcmp(detected,'-') && strcmp(overlap_char,'-')
                    Final_eqn = [Final_eqn '='];
                    
                    prev_centroid_y_coord = 1/2*(centroids(i,2)+centroids(i+1,2));
                end
                i = i+1; 
            end
    else
            
            overlap_indexs(i) = true;
            overlap_indexs;
            overlap_chars = chars(overlap_indexs);
            overlap_centroids = centroids(overlap_indexs,:)';
            overlap_boxes = boundingboxes(overlap_indexs,:)';
            bar_idx = [];
            bar_width = [];
            limit_idx = [];
            limit_height = [];
            for overlap_i = 1:length(overlap_chars)
                overlap_i;
                if strcmp(overlap_chars(overlap_i).char,'-')
                    bar_idx(end+1) = overlap_i;
                    bar_width(end+1) = overlap_boxes(3,overlap_i);
                elseif any(strcmp(overlap_chars(overlap_i).char,limit_control))
                    limit_idx(end+1) = overlap_i;
                    limit_height(end+1) = overlap_boxes(4,overlap_i);
                end
            end
            
            topleft_list = overlap_boxes(1,:);
            topright_list = overlap_boxes(1,:)+overlap_boxes(3,:);
            total_width = max(topright_list)-min(topleft_list);
            
            
            bar_idx
            bar_width
            total_width
            bar_idx = bar_idx(bar_width >= 0.8*total_width)
            
            
            if length(bar_idx) == 1
                prev_fraction = true;
                frac_bar = overlap_chars(bar_idx);
                frac_bar_centroid=overlap_centroids(:,bar_idx);
                frac_bar_box=overlap_boxes(:,bar_idx)';
                frac_y_coord = frac_bar_centroid(2);
                frac_bar_height = frac_bar_box(4);
                
                
                if ~isnan(prev_centroid_y_coord) && frac_y_coord < prev_centroid_y_coord - 7*frac_bar_height
                    is_super = true;
                    Final_eqn = [strtrim(Final_eqn) '^{'];
                end
                
                prev_centroid_y_coord = frac_y_coord;
                
                
                overlap_centroids
                frac_y_coord
                numer_idx = overlap_centroids < frac_y_coord;
                numer_idx = numer_idx(2,:);
                denom_idx = overlap_centroids > frac_y_coord;
                denom_idx = denom_idx(2,:);
                
                numer_eq_struct.filename = '';
                numer_eq_struct.characters = overlap_chars(numer_idx);
                overlap_chars(numer_idx).char
                denom_eq_struct.filename = '';
                denom_eq_struct.characters = overlap_chars(denom_idx);
                blank.val=overlap_chars(numer_idx);
                overlap_boxes(:,numer_idx)
                
                numer_str = Form_equation(blank,overlap_boxes(:,numer_idx)',overlap_centroids(:,numer_idx)');
                
                
                blank.val=overlap_chars(denom_idx);
                denom_str = Form_equation(blank,overlap_boxes(:,denom_idx)',overlap_centroids(:,denom_idx)');
                
                detected = ['\frac{' numer_str '}{' denom_str '}'];
                
                if is_super
                    is_super = false;
                    detected = [detected '}'];
                end
                
                i = find(overlap_indexs,1,'last');
                
            elseif ~isempty(limit_idx)
                [~,dom_idx] = max(limit_height); 
                dom_idx = limit_idx(dom_idx);
                limit_char = overlap_chars(dom_idx);
                limit_centroid=overlap_centroids(:,dom_idx);
                limit_y_coord = limit_centroid(2);
                
                prev_centroid_y_coord = limit_y_coord;
                
               
                top_idx = overlap_centroids < limit_y_coord;
                top_idx = top_idx(2,:);
                bottom_idx = overlap_centroids > limit_y_coord;
                bottom_idx = bottom_idx(2,:);
                
                detected = ['\' limit_char.char '\limits'];
                if any(bottom_idx)
                    bottom_eq_struct.filename = '';
                    bottom_eq_struct.characters = overlap_chars(bottom_idx);
                    blank.val=overlap_chars(bottom_idx);
                    bottom_str = Form_equation(blank,overlap_boxes(:,bottom_idx)',overlap_centroids(:,bottom_idx)');
                    
                    detected = [detected '_{' bottom_str '}'];
                end
                if any(top_idx)
                    top_eq_struct.filename = '';
                    top_eq_struct.characters = overlap_chars(top_idx);
                    blank.val=overlap_chars(top_idx);
                    top_str = Form_equation(blank,overlap_boxes(:,top_idx)',overlap_centroids(:,top_idx)');
                    
                    detected = [detected '^{' top_str '}'];
                end
                
                i = find(overlap_indexs,1,'last');
            end
            
            Final_eqn = [Final_eqn detected];
    end
    i = i+1;
end

if is_super || is_sub
    Final_eqn = [strtrim(Final_eqn) '}'];
end

    function tex_char = help_get_char(detected,index, chars,dist_to_next)
        if any(strcmp(detected,control)) || any(strcmp(detected,limit_control))
            detected = ['\' detected ' '];
        end
        
        if ~isnan(dist_to_next) && help_is_letter(detected) && help_is_letter(chars(index+1).char)
            if dist_to_next >= space
                detected = [detected '\,'];
            end
        end
        tex_char = detected;
    end


end

function control = Get_symbols()
control = {
    'alpha'
    'beta'
    'gamma'
    'delta'
    'epsilon'
    'zeta'
    'eta'
    'theta'
    'iota'
    'kappa'
    'lambda'
    'mu'
    'nu'
    'xi'
    'pi'
    'rho'
    'sigma'
    'tau'
    'upsilon'
    'phi'
    'chi'
    'psi'
    'omega'

    'Alpha'
    'Beta'
    'Gamma'
    'Delta'
    'Epsilon'
    'Zeta'
    'Eta'
    'Theta'
    'Iota'
    'Kappa'
    'Lambda'
    'Mu'
    'Nu'
    'Xi'
    'Pi'
    'Rho'
    'Sigma'
    'Tau'
    'Upsilon'
    'Phi'
    'Psi'
    'Omega'
    'int'
    'rightarrow'
    'infty'
    };
end

function out =  help_is_letter(str)
out = length(str)==1 && isstrprop(str,'alpha');
end

function new_chars = Preprocess(chars,boxes)
i = 1;
new_chars(1) = chars(1);
while i <= length(chars)
    detected = chars(i).char;
    
    if strcmp(detected,'l')
        start_x_coord = chars(i).boundingbox(1);
        
        end_x_coord = chars(i).boundingbox(3)*3 + start_x_coord;
        following_idx = find(boxes(1,:)>start_x_coord & boxes(1,:) < end_x_coord);
        following = chars(following_idx);
        lim_chars_idx =  (strcmp({following.char},'l')) | ...
            (strcmp({following.char},'.')) | ...
            (strcmp({following.char},'m'));
        lim_chars_idx = following_idx(lim_chars_idx);
        lim_chars(1) = chars(i);
        lim_chars(2:2+size(lim_chars_idx,2)-1) = chars(lim_chars_idx);
        lim_boxes(:,1) = boxes(:,i);
        lim_boxes(:,2:2+size(lim_chars_idx,2)-1) = boxes(:,lim_chars_idx);
        if length(lim_chars) == 4 % 'l','.','m'
            new_ul_x_coord = min(lim_boxes(1,:));
            new_ur_x_coord = max(lim_boxes(1,:) + lim_boxes(3,:));
            new_width = new_ur_x_coord-new_ul_x_coord;
            new_ul_y_coord = min(lim_boxes(2,:));
            new_bottom_y_coord = max(lim_boxes(2,:) + lim_boxes(4,:));
            new_height = new_bottom_y_coord - new_ul_y_coord;
            new_chars(i).char = 'lim';
            new_chars(i).boundingbox = [new_ul_x_coord new_ul_y_coord new_width new_height];
            new_chars(i).centroid = [new_ul_x_coord+new_width/2 new_ul_y_coord+new_height/2];
            
            chars(lim_chars_idx) = [];
            boxes(lim_chars_idx) = [];
        else
            
            new_chars(i) = chars(i);
        end
    else
        new_chars(i) = chars(i);
    end
    i = i+1;
    
end

end

