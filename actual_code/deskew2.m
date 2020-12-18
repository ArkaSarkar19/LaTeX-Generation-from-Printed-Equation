function [ deskew_img, deskewing_angle ] = deskew2(bw_img)
edges = edge(bw_img);
[H,T,~] = hough(edges,'RhoResolution',1,'ThetaResolution', 0.1 );

%%hough transformation
P = houghpeaks(H,4);

orientations = P(:,2);
% If all orientations are unique, take the top choice
if length(unique(orientations)) == 4    
    dominant_orientation = orientations(1);
else
    dominant_orientation = mode(P(:,2));
end

% The dominant orientation will be detected as the complement of the
% skewing angle. Subtract 90 degrees to get the deskewing angle.
% (Subtracting from 90 gets the skewing angle. Deskewing angle is the
% negative of that.)
deskewing_angle = T(dominant_orientation)-90;

% Limit angle  to -45  to 45
while(abs(deskewing_angle) > 45)
    deskewing_angle = deskewing_angle - sign(deskewing_angle)*90;
end

%% Perform the deskewing
deskew_img = imrotate(bw_img,deskewing_angle,'bilinear');

    Mrot = ~imrotate(true(size(bw_img)),deskewing_angle,'bilinear');
    deskew_img(Mrot&~imclearborder(Mrot)) = true;

end
