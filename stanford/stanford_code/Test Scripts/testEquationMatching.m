% Jim Brewer
% EE368
% Project
% Nov 30, 2015
%% Test full matching process against 8 test equation pictures. Export PDF
% of matched characters, count how many matched and weren't matched
addpath('../');
load('red_charPalette_withText_demo2.mat');
load('red_charPalette_Classifier_demo2.mat');
% Use equation_truth_mod.mat on photographs since order of detected
% elements changes.
load('equation_truth.mat');

% Get file names of images to modify to test
path = strcat(pwd,'/../Equations/Clean/');
directory = strcat(path,'*hr.jpg');
files = dir(directory);
% results=[];
% angles = -25:5:15;
% for a = 1:length(angles)
%     ang = angles(a);
    totcorrect=0;
    tot_total=0;
for i =1:length(files)
    correct = 0;
    total = 0;
    eq = imread(strcat(path,files(i).name));
    
    % Modify image "eq" for testing (rotate, scale, etc)
    % Rotate by 10deg 
    ang = 10;
    mask = true(size(eq(:,:,1)));
    mask_rot = imrotate(mask,ang);
    mask_rot = imerode(mask_rot,ones(4,4));
    eq = imrotate(eq,ang,'bilinear');
    eq(~repmat(mask_rot,1,1,3)) = 255;
    
    % Process equation to get matched characters
    eq_bin = fn_lighting_compensation(eq);
    [eq_deskew, ~] = fn_deskew2(eq_bin,true,true,5);
    eq_bin = eq_deskew;
    
    eq_chars = fn_segment(eq_bin);
    for j = 1:length(eq_chars)
       eq_chars(j).ident = fn_createIdent(eq_chars(j).img); 
    end
    % Find matches
    for j = 1:length(eq_chars)
        idx_matched = knnsearch(X_orig(:,1:length(chars(1).ident)),...
            eq_chars(j).ident,'distance','cityblock');
        eq_chars(j).char = chars(X_orig(idx_matched,end)).char;
    end
    
    % Compare found char strings vs truth char strings
    % Get truth chars
    eq_chars_truth=[];
    for k = 1:length(equation_truth)
       if(strcmp(files(i).name,equation_truth(k).filename))
           eq_chars_truth = equation_truth(k).eq_chars;
       end
    end
    
    for j = 1:min(length(eq_chars), length(eq_chars_truth))
        if(strcmp(eq_chars(j).char, eq_chars_truth(j).char))
           correct = correct + 1; 
        end
        total = total + 1;
    end
    totcorrect = totcorrect + correct;
    tot_total = tot_total + total;
    fprintf('Correct for file: %s \t %d\\%d \t %.1f%%\n',files(i).name,...
        correct,total,(correct/total)*100);
end
fprintf('\nTotal Correct: \t %d\\%d \t %.1f%%\n',...
        totcorrect,tot_total,(totcorrect/tot_total)*100);
% 
% results = [results [ang; (totcorrect/tot_total)*100]];
% fprintf('\nTotal Correct %.0f deg: \t %d\\%d \t %.1f%%\n',ang,...
%         totcorrect,tot_total,(totcorrect/tot_total)*100);
% end
% 
% figure(1);
% plot(results(1,:), results(2,:),'LineWidth',2);
% xlabel('Rotation Angle (degrees)');
% ylabel('Matched Character precentage');
% title('Matching Robustness to Equation (Clean) Rotation');
% axis([-25 15 0 100]);
% print('rotation_robustness','-dpdf');