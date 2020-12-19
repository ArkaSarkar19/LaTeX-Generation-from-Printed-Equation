function [ deskew_img, deskewing_angle ] = skew( bw_img )
%FN_DESKEW2 Deskews an image using the hough transform with optional
%background filling and edge softening.
% Optional args:
%   fill_flag: if true: pads deskew_img 
%with high values (whitespace) when
%                   rotating.
%   soften_flag: if true: Softens the image using the fn_soften_edges
%                   function.
%   soften_size: Determines the size of the operation used in
%   fn_soften_edges.


%% Create the edge map for the hough transform
BW = edge(bw_img,'canny');
figure(10);
imshow(BW);

[H,T,rho] = hough(BW);
P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));

lines = houghlines(BW,T,rho,P,'FillGap',5,'MinLength',7);

figure, imshow(bw_img), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red'); %%change color of the longest segment



orientations = P(:,2);
if length(unique(orientations)) == 4    
    dominant_orientation = orientations(1);
else
    dominant_orientation = mode(P(:,2));
end

    
deskewing_angle = T(dominant_orientation)-90;

if (deskewing_angle<0)
    deskewing_angle=mod(deskewing_angle,-90);
else
    deskewing_angle=mod(deskewing_angle,90);

    
end


deskew_img = imrotate(bw_img,deskewing_angle,'bilinear');
Mrot = ~imrotate(true(size(bw_img)),deskewing_angle,'bilinear');%to fix the background.
deskew_img(Mrot&~imclearborder(Mrot)) = true;



end