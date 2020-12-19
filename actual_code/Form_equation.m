function Final_eqn = Form_equation( characters,boundingboxes,centroids )

control = Get_symbols;
limit_control = {'int','sum','prod','lim'}; 

space = 10; 

Final_eqn = '';

chars = characters.val;
num_chars = length(chars);



boxes = zeros(4,num_chars);

for i = 1:num_chars
    boxes(:,i) = boundingboxes(i,:);
end
boxes;


[~, idxs] = sort(boxes(1,:));
chars = chars(idxs);
boxes = boxes(:,idxs);

chars = Preprocess(chars,boxes);

if length(chars) ~= num_chars
    num_chars = length(chars);
    boxes = zeros(4,num_chars);
    for i = 1:num_chars
        boxes(:,i) = chars(i).boundingbox;
    end
end

prev_centroid_y_coord = NaN;
prev_fraction = false; 
is_super = false; 
is_sub = false; 

i = 1;
while i <= num_chars

    detected = chars(i).char;
    
    %% Check for overlaps
    top_left = boxes(1,i);
    top_right = top_left+boxes(3,i);
    

    
    overlap_indexs = boxes(1,:)>=top_left & boxes(1,:)<=top_right;
    cur_overlap_boxes = boxes(:,overlap_indexs);
    cur_overlap_top = min(cur_overlap_boxes(2,:));
    cur_overlap_bottom = max(cur_overlap_boxes(2,:)+cur_overlap_boxes(4,:));
    total_height = cur_overlap_bottom - cur_overlap_top;
    overlap_indexs = boxes(1,:)>=top_left & boxes(1,:)<=top_right &...
    boxes(4,:) < 0.7*total_height;
    

    cur_overlap_boxes = boxes(:,overlap_indexs);
    topleft_list = cur_overlap_boxes(1,:);
    topright_list = cur_overlap_boxes(1,:)+cur_overlap_boxes(3,:);
    top_list = cur_overlap_boxes(2,:);
    bottom_list = cur_overlap_boxes(2,:) + cur_overlap_boxes(4,:);
    
    
    while ~isempty(topleft_list) &&...
            (min(topleft_list) < top_left || max(topright_list) > top_right)
        top_left = min(topleft_list);
        top_right = max(topright_list);
        total_height = max(bottom_list)-min(top_list);
        overlap_indexs = boxes(1,:)>=top_left & boxes(1,:)<=top_right &...
            boxes(4,:) < 0.7*total_height;
        cur_overlap_boxes = boxes(:,overlap_indexs);
        topleft_list = cur_overlap_boxes(1,:);
        topright_list = cur_overlap_boxes(1,:)+cur_overlap_boxes(3,:);
    end
    
    

    overlap_indexs(i) = false;
    num_overlaps = sum(overlap_indexs);
    
    if (num_overlaps ==0)
        
            
            if i-1 > 0 && ~prev_fraction
                ll_corner = boxes(2,i)+boxes(4,i);
                ul_corner = boxes(2,i);
                
                if ll_corner <= ceil(prev_centroid_y_coord)
                    
                    if ~strcmp(detected,'-') || ll_corner < (prev_centroid_y_coord-5*boxes(4,i))
                        
                        is_super = true;
                        Final_eqn = [strtrim(Final_eqn) '^{'];
                        
                        
                        prev_centroid_y_coord = centroids(i,2);
                        
                        if i+1<= num_chars
                            dist_to_next = boxes(1,i+1)-top_right;
                        else
                            dist_to_next = NaN;
                        end
                        detected = help_get_char(detected,i,chars,dist_to_next);
                    end
                elseif ul_corner >= floor(prev_centroid_y_coord)
                    if is_super
                        is_super = false;
                        Final_eqn = [strtrim(Final_eqn) '}'];
                        prev_fraction = false;
                        
                        prev_centroid_y_coord = centroids(i,2);
                        
                        if i+1<= num_chars
                            dist_to_next = boxes(1,i+1)-top_right;
                        else
                            dist_to_next = NaN;
                        end
                        detected = help_get_char(detected,i,chars,dist_to_next);
                    else
                        is_sub = true;
                        Final_eqn = [strtrim(Final_eqn) '_{'];
                        
                        if i+1<= num_chars
                            dist_to_next = boxes(1,i+1)-top_right;
                        else
                            dist_to_next = NaN;
                        end
                        detected = help_get_char(detected,i,chars,dist_to_next);
                        
                        prev_centroid_y_coord = centroids(i,2);
                    end
                else
                    prev_fraction = false;
                    
                    prev_centroid_y_coord = centroids(i,2);
                    
                    if i+1<= num_chars
                        dist_to_next = boxes(1,i+1)-top_right;
                    else
                        dist_to_next = NaN;
                    end
                    detected = help_get_char(detected,i,chars,dist_to_next);
                end
            else
                prev_fraction = false;
                
                prev_centroid_y_coord = centroids(i,2);
                
                if i+1<= num_chars
                    dist_to_next = boxes(1,i+1)-top_right;
                else
                    dist_to_next = NaN;
                end
                detected = help_get_char(detected,i,chars,dist_to_next);
            end
            Final_eqn = [Final_eqn detected];
    elseif (num_overlaps==1)
            assert(overlap_indexs(i+1))
            overlap_ul = boxes(1,i+1);
            if abs(overlap_ul-top_right) <= 1
                if any(strcmp(detected,control))
                    detected = ['\' detected ' '];
                end
                
                if i+1 <= num_chars && help_is_letter(detected) && help_is_letter(chars(i+1).char)
                    dist_to_next = boxes(1,i+1)-top_right;
                    if dist_to_next >= space
                        detected = [detected '\,'];
                    end
                end
                Final_eqn = [Final_eqn detected];
                
                prev_centroid_y_coord = chars(i).centroid(2);
            else
                overlap_char = chars(i+1).char;
                
                if strcmp(detected,'-') && strcmp(overlap_char,'-')
                    Final_eqn = [Final_eqn '='];
                    
                    prev_centroid_y_coord = 1/2*(centroids(i,2)+centroids(i+1,2));
                end
                i = i+1; 
            end
    else
            
            overlap_indexs(i) = true;
            overlap_indexs;
            overlap_chars = chars(overlap_indexs);
            overlap_centroids = centroids(overlap_indexs,:)';
            overlap_boxes = boundingboxes(overlap_indexs,:)';
            bar_idx = [];
            bar_width = [];
            limit_idx = [];
            limit_height = [];
            for overlap_i = 1:length(overlap_chars)
                overlap_i;
                if strcmp(overlap_chars(overlap_i).char,'-')
                    bar_idx(end+1) = overlap_i;
                    bar_width(end+1) = overlap_boxes(3,overlap_i);
                elseif any(strcmp(overlap_chars(overlap_i).char,limit_control))
                    limit_idx(end+1) = overlap_i;
                    limit_height(end+1) = overlap_boxes(4,overlap_i);
                end
            end
            
            topleft_list = overlap_boxes(1,:);
            topright_list = overlap_boxes(1,:)+overlap_boxes(3,:);
            total_width = max(topright_list)-min(topleft_list);
            
            
            bar_idx
            bar_width
            total_width
            bar_idx = bar_idx(bar_width >= 0.8*total_width)
            
            
            if length(bar_idx) == 1
                prev_fraction = true;
                frac_bar = overlap_chars(bar_idx);
                frac_bar_centroid=overlap_centroids(:,bar_idx);
                frac_bar_box=overlap_boxes(:,bar_idx)';
                frac_y_coord = frac_bar_centroid(2);
                frac_bar_height = frac_bar_box(4);
                
                
                if ~isnan(prev_centroid_y_coord) && frac_y_coord < prev_centroid_y_coord - 7*frac_bar_height
                    is_super = true;
                    Final_eqn = [strtrim(Final_eqn) '^{'];
                end
                
                prev_centroid_y_coord = frac_y_coord;
                
                
                overlap_centroids
                frac_y_coord
                numer_idx = overlap_centroids < frac_y_coord;
                numer_idx = numer_idx(2,:);
                denom_idx = overlap_centroids > frac_y_coord;
                denom_idx = denom_idx(2,:);
                
                numer_eq_struct.filename = '';
                numer_eq_struct.characters = overlap_chars(numer_idx);
                overlap_chars(numer_idx).char
                denom_eq_struct.filename = '';
                denom_eq_struct.characters = overlap_chars(denom_idx);
                blank.val=overlap_chars(numer_idx);
                overlap_boxes(:,numer_idx)
                
                numer_str = Form_equation(blank,overlap_boxes(:,numer_idx)',overlap_centroids(:,numer_idx)');
                
                
                blank.val=overlap_chars(denom_idx);
                denom_str = Form_equation(blank,overlap_boxes(:,denom_idx)',overlap_centroids(:,denom_idx)');
                
                detected = ['\frac{' numer_str '}{' denom_str '}'];
                
                if is_super
                    is_super = false;
                    detected = [detected '}'];
                end
                
                i = find(overlap_indexs,1,'last');
                
            elseif ~isempty(limit_idx)
                [~,dom_idx] = max(limit_height); 
                dom_idx = limit_idx(dom_idx);
                limit_char = overlap_chars(dom_idx);
                limit_centroid=overlap_centroids(:,dom_idx);
                limit_y_coord = limit_centroid(2);
                
                prev_centroid_y_coord = limit_y_coord;
                
               
                top_idx = overlap_centroids < limit_y_coord;
                top_idx = top_idx(2,:);
                bottom_idx = overlap_centroids > limit_y_coord;
                bottom_idx = bottom_idx(2,:);
                
                detected = ['\' limit_char.char '\limits'];
                if any(bottom_idx)
                    bottom_eq_struct.filename = '';
                    bottom_eq_struct.characters = overlap_chars(bottom_idx);
                    blank.val=overlap_chars(bottom_idx);
                    bottom_str = Form_equation(blank,overlap_boxes(:,bottom_idx)',overlap_centroids(:,bottom_idx)');
                    
                    detected = [detected '_{' bottom_str '}'];
                end
                if any(top_idx)
                    top_eq_struct.filename = '';
                    top_eq_struct.characters = overlap_chars(top_idx);
                    blank.val=overlap_chars(top_idx);
                    top_str = Form_equation(blank,overlap_boxes(:,top_idx)',overlap_centroids(:,top_idx)');
                    
                    detected = [detected '^{' top_str '}'];
                end
                
                i = find(overlap_indexs,1,'last');
            end
            
            Final_eqn = [Final_eqn detected];
    end
    i = i+1;
end

if is_super || is_sub
    Final_eqn = [strtrim(Final_eqn) '}'];
end

    function tex_char = help_get_char(detected,index, chars,dist_to_next)
        if any(strcmp(detected,control)) || any(strcmp(detected,limit_control))
            detected = ['\' detected ' '];
        end
        
        if ~isnan(dist_to_next) && help_is_letter(detected) && help_is_letter(chars(index+1).char)
            if dist_to_next >= space
                detected = [detected '\,'];
            end
        end
        tex_char = detected;
    end


end

function control = Get_symbols()
control = {
    'alpha'
    'beta'
    'gamma'
    'delta'
    'epsilon'
    'zeta'
    'eta'
    'theta'
    'iota'
    'kappa'
    'lambda'
    'mu'
    'nu'
    'xi'
    'pi'
    'rho'
    'sigma'
    'tau'
    'upsilon'
    'phi'
    'chi'
    'psi'
    'omega'

    'Alpha'
    'Beta'
    'Gamma'
    'Delta'
    'Epsilon'
    'Zeta'
    'Eta'
    'Theta'
    'Iota'
    'Kappa'
    'Lambda'
    'Mu'
    'Nu'
    'Xi'
    'Pi'
    'Rho'
    'Sigma'
    'Tau'
    'Upsilon'
    'Phi'
    'Psi'
    'Omega'
    'int'
    'rightarrow'
    'infty'
    };
end

function out =  help_is_letter(str)
out = length(str)==1 && isstrprop(str,'alpha');
end

function new_chars = Preprocess(chars,boxes)
i = 1;
new_chars(1) = chars(1);
while i <= length(chars)
    detected = chars(i).char;
    
    if strcmp(detected,'l')
        start_x_coord = chars(i).boundingbox(1);
        
        end_x_coord = chars(i).boundingbox(3)*3 + start_x_coord;
        following_idx = find(boxes(1,:)>start_x_coord & boxes(1,:) < end_x_coord);
        following = chars(following_idx);
        lim_chars_idx =  (strcmp({following.char},'l')) | ...
            (strcmp({following.char},'.')) | ...
            (strcmp({following.char},'m'));
        lim_chars_idx = following_idx(lim_chars_idx);
        lim_chars(1) = chars(i);
        lim_chars(2:2+size(lim_chars_idx,2)-1) = chars(lim_chars_idx);
        lim_boxes(:,1) = boxes(:,i);
        lim_boxes(:,2:2+size(lim_chars_idx,2)-1) = boxes(:,lim_chars_idx);
        if length(lim_chars) == 4 % 'l','.','m'
            new_ul_x_coord = min(lim_boxes(1,:));
            new_ur_x_coord = max(lim_boxes(1,:) + lim_boxes(3,:));
            new_width = new_ur_x_coord-new_ul_x_coord;
            new_ul_y_coord = min(lim_boxes(2,:));
            new_bottom_y_coord = max(lim_boxes(2,:) + lim_boxes(4,:));
            new_height = new_bottom_y_coord - new_ul_y_coord;
            new_chars(i).char = 'lim';
            new_chars(i).boundingbox = [new_ul_x_coord new_ul_y_coord new_width new_height];
            new_chars(i).centroid = [new_ul_x_coord+new_width/2 new_ul_y_coord+new_height/2];
            
            chars(lim_chars_idx) = [];
            boxes(lim_chars_idx) = [];
        else
            
            new_chars(i) = chars(i);
        end
    else
        new_chars(i) = chars(i);
    end
    i = i+1;
    
end

end

