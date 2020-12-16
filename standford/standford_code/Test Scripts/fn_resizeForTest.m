function [ ident ] = fn_resizeForTest(img, num, minScale, maxScale )
%fn_resizeForTest Creates identity vectors for the passed image for the
%   number of different scales between minScale and maxScale. ident rows
%   are each ident vector

%   TODO: Detailed explanation goes here

addpath('../');
dS = (maxScale - minScale) / (num-1);
ident = zeros(num,22); %***Width manually set from current fn_createIdent.m
for i = 1:num
    scale = dS * (i-1) + minScale;
    img_scl = imresize(img, scale);
    % Threshold to binary again
    th = graythresh(img_scl);
    img_bin = img_scl;
    img_bin(img_scl <= th) = 0;
    img_bin(img_scl > th) = 1;
    ident(i,:) = fn_createIdent(img_bin);
end

end

