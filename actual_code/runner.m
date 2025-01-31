%% 
close all;
clear
fileName = 'Clean/eq4_hr.jpg';
showFigs=true;
outputName='../test';
dir = strcat(pwd,'/Equations/');
eq = imread(strcat(dir,fileName));
figure("Name","Original Image");                                                              %%%NORMAL IMAGE
imshow(eq);
%% 
load('red_charPalette_withText_demo2.mat');
load('red_charPalette_Classifier_demo2.mat');
%% 
close all; 


bin_image=binarization(eq);
figure("Name","Binary Image");                                                              %%%binarize IMAGE
imshow(bin_image);
%% 
close all; 
[eq_deskew, ~] = skew_correction(bin_image);



%% 

close all;
[ characters_image,characters_centroids,characters_boxes ] = segmentation(eq_deskew);
characters_ident(size(characters_image,1)).ident=[];
%% 
close all;

for i = 1:length(characters_image)
   character_ident(i).ident = Find_characteristics(characters_image(i).img); 
end
%% 

close all;
for i = 1:length(characters_image)

    
    idx_matched = similarity_function(character_ident(i).ident, "Manhattan");
    character_ident(i).char = chars(X_orig(idx_matched,end)).char;
     
 
  
    
    if(true)
        figure("Name",strcat("character",int2str(i)))
        imshow(characters_image(i).img);
        title('Input');
        
       figure("Name",strcat("Matched character",int2str(i)))
        title('Matched');

        imshow(chars(X_orig(idx_matched,end)).img);
                title('Matched');

       
    end
end




%% 
close all;
blank.val=character_ident;
clc
final_string = Form_equation(blank,characters_boxes,characters_centroids);
%% 
final_string


%% 




