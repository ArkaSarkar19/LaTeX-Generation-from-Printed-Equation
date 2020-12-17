function bw_img = binarization(img)
gray_img=rgb2gray(img);
thresh = graythresh(gray_img); %uses otsu's method to give us the threshold which has the least interclass variance.
bw_img = im2bw(gray_img,thresh);
end