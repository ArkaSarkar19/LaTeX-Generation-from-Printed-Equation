% Project - Import Templates for Character Recognition and create
% Classifier matrix
% Jim Brewer
% EE368
% Nov 18, 2015

% Helper Script to create the MatLab data of the character templates and
% their character string and identifier data. Output used by "fn_main.m"
%   **Be sure to set and uncomment the names for the saved out character
%   palette and classifier matrices below.

palette = im2double(rgb2gray(imread(strcat(pwd,'/reduced_characterPalette.jpg'))));
th = graythresh(palette);
pal_bin = palette;
pal_bin(palette <= th) = 0;
pal_bin(palette > th) = 1;

%% Extract character templates
chars = fn_segment(pal_bin);

% Remove doubles from non-contigous characters (manually review and remove)
% Check for accuracy. These are for reduced_characterPalette.jpg
chars(115)=[]; %.
chars(101)=[]; %.
chars(100)=[]; %.
chars(90)=[]; %.
chars(4)=[]; %~

% Add "truth" labels (verify order in text file versus char struct)
fileID = fopen(strcat(pwd,'/reduced_characters_noslash.txt'));
text_chars = textscan(fileID,'%s', 'delimiter',',');
text_chars = char(text_chars{1,1});
fclose(fileID);

%% Create identifiers for each character
for i = 1:length(chars)
    chars(i).ident = fn_createIdent(chars(i).img);
    chars(i).char = strtrim(text_chars(i,:));
end

%% Save out character palette
% save('red_charPalette_withText.mat','chars');

%% Train Nearest Neighbor Classifier
% Create data matrix for KNN Search with just original templates
X_orig = zeros(length(chars), length(chars(1).ident) + 1);
for i = 1:length(chars)
    X_orig(i,1:length(chars(1).ident)) = chars(i).ident;
    X_orig(i,end) = i;
end
% save('red_charPalette_Classifier.mat','X_orig');