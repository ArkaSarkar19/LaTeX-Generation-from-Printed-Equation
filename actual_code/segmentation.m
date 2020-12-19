function [ characters ] = segmentation(eq, showFigs, figNum)
%fn_segment Segment and return characters from binary image



% Invert image
eq_inv = ones(size(eq)) - eq;

% Create Edge Map
se = ones(3,3);

exp = imerode(eq_inv, se);

figure, imshow(exp, []);
title("erode");

eq_edges = xor(exp, eq_inv);

figure, imshow(eq_edges, []);
title("eq edges");
% Find Centroid for x/y location of each character
s = regionprops(eq_edges, 'Centroid');
bb = regionprops(eq_edges, 'BoundingBox');

centroids = cat(1,s.Centroid);
figure, imshow(eq_edges, []);
title("centroids");
hold on
plot(centroids(:,1),centroids(:,2),'b*')
hold off


% Matrix of centroid locations
loc = cat(1, s.Centroid);
% Matrix of boundingbox corner and sizes
boundingboxes = cat(1, bb.BoundingBox);
% Round to whole pixel values and ensure still contains character for
%   character extraction purposes
boundingboxes = floor(boundingboxes);
boundingboxes(:,3:4) = boundingboxes(:,3:4) + 1;


% Init struct to contain extracted characters and their x/y locations
characters(size(loc,1)).centroid = [];


% Extract each character from binary image based on bounding boxes
%   If there are issues with overlapping characters into nearby bounding
%   boxes, can use the Convex Hulls to create masks around desired char
for i = 1 : size(loc,1)
    characters(i).centroid = loc(i,:);
    % Add Bounding Box Information to equations
    characters(i).boundingbox = boundingboxes(i,:);
    
    % Check if more than 1 object in region, if so only extract largest
    %   Assume largest is like the sqrt, etc and the others will be
    %   extracted with other centroids
    region = eq_inv(boundingboxes(i,2):boundingboxes(i,2)+boundingboxes(i,4),...
        boundingboxes(i,1):boundingboxes(i,1)+boundingboxes(i,3),:);
   
        characters(i).img = eq(boundingboxes(i,2):boundingboxes(i,2)+boundingboxes(i,4),...
        boundingboxes(i,1):boundingboxes(i,1)+boundingboxes(i,3),:);
    
    
end
dummy=zeros(size(loc,1));
for i =1:size(dummy)
    for j=1:size(dummy)
        if (i~=j)
            if ((boundingboxes(i,2)<=boundingboxes(j,2))&&(boundingboxes(i,2)+boundingboxes(i,4)>boundingboxes(j,2)+boundingboxes(j,4))&&(boundingboxes(i,1)<=boundingboxes(j,1))&&(boundingboxes(i,1)+boundingboxes(i,3)>boundingboxes(j,1)+boundingboxes(j,3)))
                dummy(j)=1
            end

        end
    end
    

end
for i=1:size(dummy)
    if (dummy(i)==1)
        characters(i)=[];
        boundingboxes(i,:)=[];
        loc(i,:)=[];
    end
    
        
end



%% Only show figures if boolean passed to show them
if showFigs
    % Plot centroids, bounding boxes, convex hulls on image
    figure(figNum);
    imshow(eq_edges);
    for i = 1:size(loc, 1) % Assumes same # centroid and convex hull
        % Plot centroids
        x = loc(i,1);
        y = loc(i,2);
        text(x, y, '*' ,'Color', 'yellow', 'FontSize', 14);
        
        % Plot Bounding Box
        rectangle('position',boundingboxes(i,:),'Edgecolor','g')
    end

    % Show all extracted characters
    figure(figNum + 1);
    for i =1:size(characters,2)
        dim = ceil(sqrt(size(characters,2)));
       subplot(dim, dim,i);
       imshow(characters(i).img);
       title(sprintf('%d',i));
    end
end

end