function [ characters_image,characters_centroids,characters_boxes ] = segmentation(img)





exp = imerode(ones(size(img)) - img, ones(3,3));

figure("Name","Eroded image"), imshow(exp, []);

img_edges = xor(imerode(ones(size(img)) - img, ones(3,3)), ones(size(img)) - img);



Ctr = regionprops(img_edges, 'Centroid');
BoundingB = regionprops(img_edges, 'BoundingBox');



centroids = cat(1, Ctr.Centroid);

boundingboxes = cat(1, BoundingB.BoundingBox);
boundingboxes = floor(boundingboxes);
boundingboxes(:,3:4) = boundingboxes(:,3:4) + 1;


characters_image(size(centroids,1)).img=[];

dummy=zeros(size(centroids,1),1);
for i =1:size(dummy,1)
    for j=1:size(dummy,1)
        if (i~=j)
            if ((boundingboxes(i,2)<=boundingboxes(j,2))&&(boundingboxes(i,2)+boundingboxes(i,4)>boundingboxes(j,2)+boundingboxes(j,4))&&(boundingboxes(i,1)<=boundingboxes(j,1))&&(boundingboxes(i,1)+boundingboxes(i,3)>boundingboxes(j,1)+boundingboxes(j,3)))
                dummy(j)=1;
            end

        end
    end
    

end
counter=0;
for i=1:size(dummy)
    if (dummy(i)==1)
        characters_image(i-counter)=[];
        boundingboxes(i-counter,:)=[];
        centroids(i-counter,:)=[];
        counter=counter+1;
    end
    
        
end



for i = 1 : size(centroids,1)
    
    characters_image(i).img = img(boundingboxes(i,2):boundingboxes(i,2)+boundingboxes(i,4),boundingboxes(i,1):boundingboxes(i,1)+boundingboxes(i,3),:);
    
    
end



figure("Name","Centroids");

imshow(img_edges);

for i = 1:size(centroids, 1) 
    hold on ;
    x = centroids(i,1);
    y = centroids(i,2);
    text(x, y, '*' ,'Color', 'red', 'FontSize', 10);
    hold off;

end
figure("Name","Bounding Boxes");

imshow(img_edges);

for i = 1:size(boundingboxes, 1) 
    hold on ;
    rectangle('position',boundingboxes(i,:),'Edgecolor','y')

    hold off;

end

characters_centroids=centroids;
characters_boxes=boundingboxes;
end