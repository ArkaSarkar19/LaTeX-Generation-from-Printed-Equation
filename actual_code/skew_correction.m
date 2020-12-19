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