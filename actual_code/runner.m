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
