function [ identifier ] = create_indent(char, printIdent, showFigs)
%createIdent - Create the identifier vector with normalized central moment
% of inertia, the number of chracter crossings for each concentric circle,
% and the 2 largest arc length ratio for k-1 largest concentric circles.
% Also uses the other 6 Hu moments as part of the identifier.
%   
%   char - Binary image of character to create identifier for. Assumes
%       character is black.
%   printIdent - Prints out central moment, counts and ratio identifier
%   showFigs - Set to true to show extracted circles and figure with
%       circles overlaid.
%
%   "coding" is a struct that contains the extracted vectors, counts, and 
%       and rations but is not currently returned.
%   
%   Assumes: Character is binary image (character = 0) tightly cropped.

if nargin == 1
    printIdent = false;
    showFigs = false;
elseif nargin == 2
    showFigs = false;
end

% Parameters
k = 8;
% IN2, R1-R7, D2-D8, 6 Hu Moments
identifier = zeros(1,2*k+6);

char_inv = ones(size(char)) - char;

figure, imshow(char, []);
title("char")

figure, imshow(char_inv, []);
title("char_inv")
% If white border on any edge, pad with more 1s
pad = 3;
% top
if(sum(char_inv(1,:)) == 0)
    char = [ones(pad,size(char,2)) ;char];
end
% bottom
if(sum(char_inv(end,:)) == 0)
    char = [char; ones(pad,size(char,2))];
end
% left
if(sum(char_inv(:,1)) == 0)
    char = [ones(size(char,1),pad) char];
end
% right
if(sum(char_inv(:,end)) == 0)
    char = [char ones(size(char,1),pad)];
end
% Recalculate for centroid
char_inv = ones(size(char)) - char;

% Total number of pixels in character
N = sum(~char(:));
% Find x,y positions of 0 elements (character)
[y, x] = find(~char);

% Find centroid
cent = regionprops(char_inv,'centroid');
cent = cat(1, cent.Centroid);
cent_x = cent(1);
cent_y = cent(2);
cent_x
% Find central moment of inertia
I = sum((x - cent_x).^2 + (y - cent_y).^2);
IN2 = I / N^2;
identifier(1) = IN2;
if(printIdent)
   fprintf('Normalized Central Moment: %.4f\n', IN2); 
end

% Extract Circular Topology Template and Diameter Ratio
y = size(char,1);
x = size(char,2);
dr = max(max(x - cent_x,cent_x), max(y - cent_y,cent_y)) / (k+1);
[c, r] = meshgrid(1:x, 1:y);

% Display Extracted circles
if(showFigs)
    figure(1);
    % Create copy of character to display circles overlaid
    tempcirc = char;
end

% Init coding vector for storing extra information
coding(k).circVec=[];
for i = 1:k
    rad = dr * i;
    % Create circle line of logical 1s to extract template
    C = xor(sqrt((c-cent_x).^2+(r-cent_y).^2)<=rad, ...
        sqrt((c-cent_x).^2+(r-cent_y).^2)<=(rad-1));
    % Extract circular template linear indices (sort based on theta)
    cidx = find(C);
    % Create matrix linear indices, corresponding x and y values
    vals = [cidx c(cidx)-cent_x r(cidx)-cent_y zeros(size(cidx,1),1)];
    % Add corresponding theta values
    [TH, ~] = cart2pol(vals(:,2), vals(:,3));
    vals(:,4) = TH;
    % Sort based on theta values so in order around circle
    [~, order] = sort(vals(:,4));
    sortedvals = vals(order,:);
    % Use linear indices to extract desired values in correct theta order
    circVec = char(sortedvals(:,1));
    
    % Show circles overlaid on character
    if(showFigs)
        tempcirc(C) = .5;
    end
    
    % Extract values from character based on circular indices. Start at
    %   0deg going CCW.
    circVec = imopen(circVec,ones(2,1));
    
    % Save vector for debugging or future use
    coding(i).circ = circVec;
    
    % If character is so small nothing is extracted, leave ident as 0
    if(~isempty(circVec))
        % Count number of sections of at least 2 0s. This is the number of
        % character sections the circle goes through.
        cnt = strfind([1 1 circVec'],[0 0]);
        coding(i).count = length(cnt(diff([1 cnt])~=1));

        % Save count to identifier vector
        identifier(1+i) = coding(i).count;
        if(showFigs)
            subplot(2,4,i);
            imshow(coding(i).circ);
            str = sprintf('%d',coding(i).count);
            title(str);
        end

        % Find 2 longest arcs of background and divide diff by total length
        % If 1s (background) wrap around, move all 1s to one side.
        circ = length(circVec);
        if(circVec(1) == 1 && circVec(end) == 1)
            idx = find(circVec==0,1,'last');
            circVec = [circVec(idx+1:end); circVec(1:idx)];
        end
        % Vector of length of each section of 1s
        B = [0 circVec' 0]; % Pad with 0s for diff
        bgrd_len = find(diff(B)==-1) - find(diff(B)==1);
        % If no crossings, set ratio to 0
        if(coding(i).count < 1);
            d2 = 0;
            d1 = 0;
        else
            % Get two longest sections (arcs)
            [d2, idx_d2] = max(bgrd_len);
            bgrd_len(idx_d2)=NaN;
            d1 = max(bgrd_len);
        end
        % Find ratio of difference of two longest arcs by circumference
        coding(i).ratio = (d2-d1) / circ;
        if(i > 1) % Only keep k-1 largest ratios
            if(~isnan(coding(i).ratio))
                identifier(k+i) = coding(i).ratio;
            end
        end
        if(printIdent)
            fprintf('Count=%d, d2=%d, d1=%d, c=%d, ratio = %.4f\n',...
                coding(i).count,d2,d1,circ,coding(i).ratio);
        end
    end
end

% Add Hu Moments to end of Identifier (work on inverse so char is 0)
eta_mat = SI_Moment(char_inv) ;
hu_arr = Hu_Moments(eta_mat);
% 1st moment is IN2, above, in first slot
identifier(2*k+1:end) = hu_arr(2:end);

% Show character with circles overlaid
if(showFigs)
    figure(2);
    imshow(tempcirc);
end

end
