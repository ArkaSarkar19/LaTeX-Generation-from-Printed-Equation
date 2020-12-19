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