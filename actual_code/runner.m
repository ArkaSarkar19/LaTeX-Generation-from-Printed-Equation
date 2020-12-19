%% 
fileName = 'Clean/R2.jpg';
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
sum(sum(eq_bin))
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

for i = 1:length(eq_chars)
   eq_chars(i).ident = create_indent(eq_chars(i).img); 
end

%% 
for i = 1:length(eq_chars)
    idx_matched = knnsearch(X_orig(:,1:length(chars(1).ident)),...
            eq_chars(i).ident,'distance','cityblock');
    % Set the matched character value
    eq_chars(i).char = chars(X_orig(idx_matched,end)).char;
    
    % Code to show input characters and their determined matches
    if(showFigs)
        figure(6);
        subplot(2,length(eq_chars),i);
        imshow(eq_chars(i).img);
        title('Input');
        subplot(2,length(eq_chars),i+length(eq_chars));
        imshow(chars(X_orig(idx_matched,end)).img);
        if(eq_chars(i).char(1) == '\')
            printChar = strcat('\',eq_chars(i).char);
        else
           printChar = eq_chars(i).char;
        end
        str = sprintf('Match: %s',printChar);
        title(str);
    end
end

%% Pass struct of segmented characters (eq_chars) with matched character 
% data to equation creator

EqStruct.characters = eq_chars;
eq_string = fn_assemble_eq(EqStruct);

%% Output LaTeX Code
writeTex(eq_string, strcat(dir,outputName));
