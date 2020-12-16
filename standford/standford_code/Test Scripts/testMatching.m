% Test Classifier/Matching
% EE368
% Nov 20, 2015

% Assume character palette has been created and saved out. If not, use
% importCharacterTemplates to do so
addpath('../');
load('red_charPalette_withText_demo2.mat');

%% Create identifiers for each character
for i = 1:length(chars)
    chars(i).ident = fn_createIdent(chars(i).img);
    % 0 out count for how many times that character is a false positive
    chars(i).wrong = 0;
end

%% Train Nearest Neighbor Classifier
% Create data matrix for KNN Search with just original templates
X_orig = zeros(length(chars), length(chars(1).ident) + 1);
for i = 1:length(chars)
    X_orig(i,1:length(chars(1).ident)) = chars(i).ident;
    X_orig(i,end) = i;
end
% save('red_charPalette_Classifier.mat','X_orig');

%% Find percentage correct with resized/rotated test images over entire palette
correct = 0;
total = 0;
% Scale Robustness Test
scale = true;
numResize = 4;
maxScale = 1.25;
minScale = .5;
% Skew Robustness Test
rotate = false;
numRot = 7;
minAng = -15;
maxAng = 15;
scaleMatches = zeros(1,numResize);
rotMatches = zeros(1,numRot);
for j = 1:length(chars)
    if(scale)
        ident_testing = fn_resizeForTest(chars(j).img, numResize, ...
            minScale, maxScale);
    elseif(rotate)
        ident_testing = fn_rotateForTest(chars(j).img, numRot, ...
            minAng, maxAng);
    end
    chars(j).match = [];
    for i = 1:size(ident_testing,1)
       idx_hat = knnsearch(X_orig(:,1:length(chars(1).ident)),...
           ident_testing(i,:),'distance','cityblock');
       % Record what it thinks the match is
       chars(j).match = [chars(j).match X_orig(idx_hat,end)];
       if(j == X_orig(idx_hat,end))
           correct = correct + 1;
           % Increment the correct number for that scale
           if(scale)
               scaleMatches(i) = scaleMatches(i) + 1;
           elseif(rotate)
               rotMatches(i) = rotMatches(i) + 1;
           end
       else
          % Mark which ones are common false matches
          chars(j).wrong = chars(j).wrong + 1;
       end
       total = total + 1;
    end
end

if(scale)
    fprintf('Resized Test Images (%d) with KNNSEARCH:\n',numResize);
    fprintf('# of Templates:%d\n',length(chars));
    fprintf('# of Test Images:%d\n\n',total);
    for i = 1:numResize
       fprintf('Correct for Scale %.2f: %.0f \t %.1f%%\n',...
            ((maxScale - minScale) / (numResize-1)) * (i-1) + minScale,...
            scaleMatches(i), (scaleMatches(i)/length(chars))*100);
    end
    fprintf('\nTotal Correct:\t\t%.0f \t %.1f%%\n',correct,(correct/total)*100);
elseif(rotate)
    fprintf('Rotated Test Images (%d) (CCW) with KNNSEARCH:\n',numRot);
    fprintf('# of Templates:%d\n',length(chars));
    fprintf('# of Test Images:%d\n\n',total);
    for i = 1:numRot
       fprintf('Correct for Rotation Angle %.0fDeg: %.0f \t %.1f%%\n',...
            ((maxAng - minAng) / (numRot - 1)) * (i-1) + minAng,...
            rotMatches(i), (rotMatches(i)/length(chars))*100);
    end
    fprintf('\nTotal Correct:\t\t%.0f \t %.1f%%\n',correct,(correct/total)*100);
end

%% Show which characters are being mismatched with which

matches = [(1:length(chars))' cat(1,chars.match)];
matches_bool = matches(:,2:end) ~= repmat(matches(:,1),1,numRot);
scl = 2; %Pick a scale number (1 to numResize above)
[row, col] = find(matches_bool(:,scl));
figure(1);
for i = 1:length(row)
   subplot(2,length(row),i);
   imshow(chars(row(i)).img);
   subplot(2,length(row),i+length(row));
   imshow(chars(matches(row(i),scl+1)).img);
end

%% Test with input equation
directory = strcat(pwd,'/../Equations');
% Test with simple binarization
% directory = strcat(pwd,'/Lighting');
% eq = im2double(rgb2gray(imread(strcat(directory,'/eq1_hr.jpg'))));
% th = graythresh(eq);
% eq_bin = eq;
% eq_bin(eq <= th) = 0;
% eq_bin(eq > th) = 1;

eq = imread(strcat(directory,'/Clean/eq1_hr.jpg'));

% Rotate clean images for testing
% ang = 20;
% mask = true(size(eq(:,:,1)));
% mask_rot = imrotate(mask,ang);
% mask_rot = imerode(mask_rot,ones(4,4));
% eq = imrotate(eq,ang,'bilinear');
% eq(~repmat(mask_rot,1,1,3)) = 255;

eq_bin = fn_lighting_compensation(eq);
[eq_deskew, ~] = fn_deskew2(eq_bin,true,true,3);
eq_bin = eq_deskew;
figure(1);
imshow(eq_bin);

eq_chars = fn_segment(eq_bin);
for i = 1:length(eq_chars)
   eq_chars(i).ident = fn_createIdent(eq_chars(i).img); 
end

figure(2);
for i = 1:length(eq_chars)
    subplot(2,length(eq_chars),i);
    imshow(eq_chars(i).img);
    title('Input');
    
    idx_matched = knnsearch(X_orig(:,1:length(chars(1).ident)),...
        eq_chars(i).ident,'distance','cityblock');
    eq_chars(i).char = chars(X_orig(idx_matched,end)).char;
    subplot(2,length(eq_chars),i+length(eq_chars));
    imshow(chars(X_orig(idx_matched,end)).img);
    str = sprintf('Matched %d',X_orig(idx_matched,end));
    title(str);
end



