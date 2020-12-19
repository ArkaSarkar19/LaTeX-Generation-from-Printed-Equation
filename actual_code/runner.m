%% 
close all;
fileName = 'Clean/eq5_hr.jpg';
showFigs=true;
outputName='../test';
dir = strcat(pwd,'/Equations/');
eq = imread(strcat(dir,fileName));
figure(1);                                                              %%%NORMAL IMAGE
imshow(eq);
%% 
load('red_charPalette_withText_demo2.mat');
load('red_charPalette_Classifier_demo2.mat');
%% 

eq_bin=binarization(eq);
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
[ characters_image,characters_centroids,characters_boxes ] = segmentation(eq_deskew,true,4);
characters_ident(size(characters_image,1)).ident=[];

for i = 1:length(characters_image)
   character_ident(i).ident = create_indent(characters_image(i).img); 
end

%% 
% close all;
% for i = 1:length(characters_image)
%     idx_matched = knnsearch(X_orig(:,1:length(character_ident(1).ident)),...
%             character_ident(i).ident,'distance','cityblock');
%     % Set the matched character value
%     character_ident(i).char = chars(X_orig(idx_matched,end)).char;
%      
%  
%   
%     % Code to show input characters and their determined matches
%     if(true)
%         figure(6);
%         subplot(2,length(characters_image),i);
%         imshow(characters_image(i).img);
%         title('Input');
%         subplot(2,length(characters_image),i+length(characters_image));
%         imshow(chars(X_orig(idx_matched,end)).img);
%         if(character_ident(i).char(1) == '\')
%             printChar = strcat('\',character_ident(i).char);
%         else
%            printChar = character_ident(i).char;
%         end
%         str = sprintf('Match: %s',printChar);
%         title(str);
%     end
% end

%% Pass struct of segmented characters (eq_chars) with matched character 
% data to equation creator
% 
% EqStruct.characters = eq_chars;
% eq_string = fn_assemble_eq(EqStruct);

%% Output LaTeX Code
% writeTex(eq_string, strcat(dir,outputName));
% clc

%% 


close all;
for i = 1:length(characters_image)
%     idx_matched = knnsearch(X_orig(:,1:length(character_ident(1).ident)),...
%             character_ident(i).ident,'distance','cityblock');
    % Set the matched character value
    
    idx_matched = similarity_function(character_ident(i).ident, "Manhattan");
    character_ident(i).char = chars(X_orig(idx_matched,end)).char;
     
 
  
    % Code to show input characters and their determined matches
    if(true)
        figure(6);
        subplot(2,length(characters_image),i);
        imshow(characters_image(i).img);
        title('Input');
        subplot(2,length(characters_image),i+length(characters_image));
        imshow(chars(X_orig(idx_matched,end)).img);
        if(character_ident(i).char(1) == '\')
            printChar = strcat('\',character_ident(i).char);
        else
           printChar = character_ident(i).char;
        end
        str = sprintf('Match: %s',printChar);
        title(str);
    end
end



% similarity_function(characters_image, "Manhattan");


