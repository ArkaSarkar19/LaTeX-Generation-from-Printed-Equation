%% 
fileName = 'Clean/eq1_hr.jpg';
showFigs=true;
outputName='../test';
dir = strcat(pwd,'/Equations/');
eq = imread(strcat(dir,fileName));
figure(1);                                                              %%%NORMAL IMAGE
imshow(eq);
%% 

eq_bin = binarization(eq);
figure(2);                                                              %%%binarize IMAGE
imshow(eq_bin);
%% 
close all; 
[eq_deskew, ~] = skew(eq_bin);

if(showFigs)
    figure(3);
    imshow(eq_deskew);                                                      %%%Optimize page 
end

%% Segment Equation Characters and Create Identifier

close all;
eq_chars = segmentation(eq_deskew,true,4);

% for i = 1:length(eq_chars)
%    eq_chars(i).ident = fn_createIdent(eq_chars(i).img); 
% end