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
edges = edge(bw_img);
[H,T,~] = hough(edges,'RhoResolution',1,'ThetaResolution', 0.1 );


%% Current Implementation: Taking the orientation that appears most in the 
%                           top four peaks
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

%% Perform background filling if needed
Mrot = ~imrotate(true(size(bw_img)),deskewing_angle,'bilinear');
deskew_img(Mrot&~imclearborder(Mrot)) = true;


%% Perform edge softening 
% Only if requested and only if the deskewing angle is large.
if abs(deskewing_angle) > 5
    deskew_img = fn_soften_edges(deskew_img, 5);
end


end