function bw_img = fn_lighting_otsu(img)
%FN_LIGHTING_COMPENSATION. ONLY FOR TESTING.
% Takes in an rgb image of an equation and returns a binarized version of
% the image for which uneven lighting has been compensated.
% Output will not be inverted.
addpath('../');
% Check if image needs lighting compensation by checking to see what
% proportion of pixels are in mid grayscale values
gray_img=rgb2gray(img);
[height, width] = size(img(:,:,1));
num_mid_gray = sum(sum(gray_img<240 & gray_img >15));
if num_mid_gray < 0.1*height*width % Should just be a scanned image
    thresh = graythresh(gray_img);
    bw_img = im2bw(gray_img,thresh);
else % Perform lighting compensation
    
    if height*width < 2000*1000
        blur_flag = false;
    else
        blur_flag = true;
    end
    
    if blur_flag
        fg = fspecial('gaussian', 10, 3);
        gray_img = imfilter(gray_img, fg, 'conv','replicate');
    end
    
    % Set the window size of the filter based on image dimensions
    win_size = round(min(height/60,width/60));
    [h,w] = size(gray_img);
    bw_img = false(size(gray_img));
    gray_img = im2double(gray_img);
    for i = 1+floor(win_size/2):win_size:h-ceil(win_size/2)
        for j = 1+floor(win_size/2):win_size:w-ceil(win_size/2)
            window =  gray_img(max(1,i-win_size):min(h,i+win_size),...
                                max(1,j-win_size):min(w,j+win_size));
            var_win = var(window(:));
            if var_win > 0.001
                thresh = graythresh(window);
                window = im2bw(window, thresh);
                bw_img(max(1,i-win_size):min(h,i+win_size),...
                                max(1,j-win_size):min(w,j+win_size)) = window;
            else
                mean_val = mean(window(:));
                bw_img(max(1,i-win_size):min(h,i+win_size),...
                                max(1,j-win_size):min(w,j+win_size)) = mean_val > 0.5;
            end
        end
    end
    
    % Remove small noise pixels.
    noise_size = round(0.0001*height*width);
    bw_img = bwareaopen(bw_img, noise_size);
    
    % Close gaps in edges
    se = strel('square',4);
    bw_img = imclose(bw_img,se);
    
    
    % Fill small holes (less than 5% of area of image)
    small_hole_thresh = round(0.0001*height*width);
    filled = imfill(bw_img,'holes');
    holes = filled & ~bw_img;
    lg_holes = bwareaopen(holes,small_hole_thresh);
    sm_holes = holes &~lg_holes;
    bw_img = bw_img | sm_holes;
    
    % Return image to original polarity.
    bw_img = ~bw_img;
end
end