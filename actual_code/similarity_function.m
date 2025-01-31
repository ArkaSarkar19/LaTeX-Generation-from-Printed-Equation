function character = similarity_function(img_vector, type)

if(type == "Manhattan")
    
    Orig_vectors = load('red_charPalette_Classifier_demo2.mat').X_orig;
    character_labels = load('red_charPalette_withText_demo2.mat');
    
    [r,c] = size(Orig_vectors);
    
    min_distance = 9999;
    min_label = 0;
    
    for i=1:r
        arr = Orig_vectors(i,1:end-1);
        label = Orig_vectors(i,end);
        dist = sum(abs(arr - img_vector));
        if(dist < min_distance ) 
            min_distance = dist;
            min_label = label;
        end
    end
end


if(type == "Eucledian")
    
    Orig_vectors = load('red_charPalette_Classifier_demo2.mat').X_orig;
    character_labels = load('red_charPalette_withText_demo2.mat');
    
    [r,c] = size(Orig_vectors);
    
    min_distance = 9999;
    min_label = 0;
    
    for i=1:r
        arr = Orig_vectors(i,1:end-1);
        label = Orig_vectors(i,end);
        dist = sum((arr - img_vector).^2);
        if(dist < min_distance ) 
            min_distance = dist;
            min_label = label;
        end
    end
end

character = min_label;
end