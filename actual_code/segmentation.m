function [ characters_image,characters_centroids,characters_boxes ] = segmentation(eq, showFigs, figNum)
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



% Matrix of centroid locations
centroids = cat(1, s.Centroid)
% Matrix of boundingbox corner and sizes
boundingboxes = cat(1, bb.BoundingBox);
boundingboxes = floor(boundingboxes);
boundingboxes(:,3:4) = boundingboxes(:,3:4) + 1;


% Init struct to contain extracted characters and their x/y locations
characters_image(size(centroids,1)).img=[];

dummy=zeros(size(centroids,1),1)
for i =1:size(dummy,1)
    for j=1:size(dummy,1)
        if (i~=j)
            if ((boundingboxes(i,2)<=boundingboxes(j,2))&&(boundingboxes(i,2)+boundingboxes(i,4)>boundingboxes(j,2)+boundingboxes(j,4))&&(boundingboxes(i,1)<=boundingboxes(j,1))&&(boundingboxes(i,1)+boundingboxes(i,3)>boundingboxes(j,1)+boundingboxes(j,3)))
                dummy(j)=1
            end

        end
    end
    

end
counter=0
for i=1:size(dummy)
    if (dummy(i)==1)
        characters_image(i-counter)=[];
        boundingboxes(i-counter,:)=[];
        centroids(i-counter,:)=[];
        counter=counter+1;
    end
    
        
end



for i = 1 : size(centroids,1)
    
    characters_image(i).img = eq(boundingboxes(i,2):boundingboxes(i,2)+boundingboxes(i,4),boundingboxes(i,1):boundingboxes(i,1)+boundingboxes(i,3),:);
    
    
end



%% Only show figures if boolean passed to show them
if showFigs
    for i = 1:size(centroids, 1) 
        figure();
        imshow(eq_edges);
        
        rectangle('position',boundingboxes(i,:),'Edgecolor','g')
        x = centroids(i,1);
        y = centroids(i,2);
        text(x, y, '*' ,'Color', 'yellow', 'FontSize', 14);
        
        
    end

end
characters_centroids=centroids;
characters_boxes=boundingboxes;
end