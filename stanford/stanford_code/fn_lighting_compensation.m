function bw_img = fn_lighting_compensation(img)
%FN_LIGHTING_COMPENSATION Takes in an rgb image of an equation and
% returns a binarized version of the image for which uneven lighting has
% been compensated. Output will not be inverted.

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
    
    % Get the window mean of each pixel by filtering using an averaging filter
    window_means=imfilter(gray_img,fspecial('average',win_size),'replicate');
    
    % Remove the mean and threshold. Also inverts the image.
    demeaned=window_means-gray_img-10;
    bw_img=im2bw(demeaned,0);
    
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