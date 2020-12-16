function [ ident ] = fn_rotateForTest(img, num, minAng, maxAng)
%fn_rotateForTest Rotates the passed image from minAnd to MaxAng using the
%number of increments given by "num" and passes back the indentifier vector
%for each image as a row in matrix containing them all. Positive values
%rotate CCW

%   TODO: Detailed explanation goes here
addpath('../');
dR = (maxAng - minAng) / (num-1);
ident = zeros(num,22); %***Width manually set from current fn_createIdent.m
for i = 1:num
     ang = dR * (i-1) + minAng;
     % Rotate inverse so added portion from rotation can be inverted back
     % to white
    img_rot = imrotate(ones(size(img))-img, ang);
    % Don't need ot threshold if using default "nearest neighbor" method
    ident(i,:) = fn_createIdent(ones(size(img_rot))-img_rot);
end

end

