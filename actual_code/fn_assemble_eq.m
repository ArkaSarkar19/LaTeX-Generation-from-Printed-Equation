function eq_string = fn_assemble_eq( EqStruct )
%AssembleEquation Assembles a string representation of an equation given a
% struct of extracted character information. The output eq_string will be
% in LaTeX form.
%   Input: EquationStruct is a struct returned from the fn_segment
%   function. The struct will conatain 2 fields:
%       filename    - The name of the image with the equation. Not currently
%                       used in this algorithm
%       characters  - A struct array, with each struct representing a
%                       segmented and recognized character.
%   Within the structs in the characters struct array are the following
%   fileds:
%       centroid    - The centroid of the character
%       boundingbox - the bounding box of the character
%       img         - The image of the segmented character region. Not
%                       currently used in this algorithm
%       char        - The recognized LaTeX character.

% Create control sequence recognition tables
control = help_create_control;
limit_control = {'int','sum','prod','lim'}; % Controls that can have the \limits modifier

% Set white space threshold in pixels for letters to know when to enter a
% thin space.
space = 10;

% Initialize the eq_string
eq_string = '';

% Extract the character array.
chars = EqStruct.characters;
num_chars = length(chars);


% Get boundingbox info
% Bounding boxes format:
% Upper left x, Upper left y, width, height
boxes = zeros(4,num_chars);
for i = 1:num_chars
    boxes(:,i) = chars(i).boundingbox;
end

% Make sure characters are sorted by upper left bound
[~, idxs] = sort(boxes(1,:));
chars = chars(idxs);
boxes = boxes(:,idxs);

% Do preprocessing on chars
chars = help_preprocess_chars(chars,boxes);
if length(chars) ~= num_chars
    num_chars = length(chars);
    boxes = zeros(4,num_chars);
    for i = 1:num_chars
        boxes(:,i) = chars(i).boundingbox;
    end
end

% variables for super/subscript checks
prev_centroid_y_coord = NaN;
prev_fraction = false; % Don't check for exponents if previously saw a fraction. Fractional exponents will be handled recursively
is_super = false; % Flag to know when we are in an exponent
is_sub = false; % Flag for subscript

% Loop through chars and create string
i = 1;
while i <= num_chars
    detected = chars(i).char;
    
    %% Check for overlaps
    ul_x_coord = boxes(1,i);
    ur_x_coord = ul_x_coord+boxes(3,i);
    
    % Get overlaps by checking 2 things:
    %   1. If upper left x coord falls within bounding box of overlap
    %   2. If height is significantly smaller than total height of overlap
    %   ('dt' is the main reason for check #2)
    overlap_idx = boxes(1,:)>=ul_x_coord & boxes(1,:)<=ur_x_coord;
    cur_overlap_boxes = boxes(:,overlap_idx);
    cur_overlap_top = min(cur_overlap_boxes(2,:));
    cur_overlap_bottom = max(cur_overlap_boxes(2,:)+cur_overlap_boxes(4,:));
    total_height = cur_overlap_bottom - cur_overlap_top;
    overlap_idx = boxes(1,:)>=ul_x_coord & boxes(1,:)<=ur_x_coord &...
        boxes(4,:) < 0.7*total_height;
    
    % keep checking for changed ul and ur coordinates to get all overlaps
    cur_overlap_boxes = boxes(:,overlap_idx);
    ul_list = cur_overlap_boxes(1,:);
    ur_list = cur_overlap_boxes(1,:)+cur_overlap_boxes(3,:);
    top_list = cur_overlap_boxes(2,:);
    bottom_list = cur_overlap_boxes(2,:) + cur_overlap_boxes(4,:);
    while ~isempty(ul_list) &&...
            (min(ul_list) < ul_x_coord || max(ur_list) > ur_x_coord)
        ul_x_coord = min(ul_list);
        ur_x_coord = max(ur_list);
        total_height = max(bottom_list)-min(top_list);
        overlap_idx = boxes(1,:)>=ul_x_coord & boxes(1,:)<=ur_x_coord &...
            boxes(4,:) < 0.7*total_height;
        cur_overlap_boxes = boxes(:,overlap_idx);
        ul_list = cur_overlap_boxes(1,:);
        ur_list = cur_overlap_boxes(1,:)+cur_overlap_boxes(3,:);
    end
    
    
    % Check to see how many overlaps there are, excluding the current
    % character
    overlap_idx(i) = false;
    num_overlaps = sum(overlap_idx);
    
    %% Start the assembly
    %% special case if sqrt
    if strcmp(detected,'sqrt')
        sqrt_chars = chars(overlap_idx);
        sqrt_eq_struct.fname = '';
        sqrt_eq_struct.characters = sqrt_chars;
        sqrt_str = fn_assemble_eq(sqrt_eq_struct);
        eq_string = [eq_string '\sqrt{' sqrt_str '}']; %#ok<*AGROW>
        i = i+num_overlaps+1;
        continue
    end
    
    %% Special cases depending on number of overlaps
    switch num_overlaps
        case 0
            % Check for super and sub scripts by comparing lower left with
            % stored centroid_height, unless previous character was a fraction
            if i-1 > 0 && ~prev_fraction
                ll_corner = boxes(2,i)+boxes(4,i);
                ul_corner = boxes(2,i);
                
                if ll_corner <= ceil(prev_centroid_y_coord)
                    % Check to see if this was a '-'. If so, it
                    %  should be pretty high in order to be an exponent
                    if ~strcmp(detected,'-') || ll_corner < (prev_centroid_y_coord-5*boxes(4,i))
                        % Remove trailing white space, append the superscript
                        % control sequence
                        is_super = true;
                        eq_string = [strtrim(eq_string) '^{'];
                        
                        % Store centroid to check for nested exponents
                        prev_centroid_y_coord = chars(i).centroid(2);
                        
                        if i+1<= num_chars
                            dist_to_next = boxes(1,i+1)-ur_x_coord;
                        else
                            dist_to_next = NaN;
                        end
                        % determine how detected is represented in LaTeX
                        detected = help_get_char(detected,i,chars,dist_to_next);
                    end
                elseif ul_corner >= floor(prev_centroid_y_coord)
                    if is_super
                        is_super = false;
                        % End of exponent
                        eq_string = [strtrim(eq_string) '}'];
                        prev_fraction = false;
                        
                        % Store centroid to check for exponents
                        prev_centroid_y_coord = chars(i).centroid(2);
                        
                        if i+1<= num_chars
                            dist_to_next = boxes(1,i+1)-ur_x_coord;
                        else
                            dist_to_next = NaN;
                        end
                        detected = help_get_char(detected,i,chars,dist_to_next);
                    else
                        % Subscript
                        is_sub = true;
                        eq_string = [strtrim(eq_string) '_{'];
                        
                        if i+1<= num_chars
                            dist_to_next = boxes(1,i+1)-ur_x_coord;
                        else
                            dist_to_next = NaN;
                        end
                        detected = help_get_char(detected,i,chars,dist_to_next);
                        
                        % Store centroid to check for exponents
                        prev_centroid_y_coord = chars(i).centroid(2);
                    end
                else
                    % Normal op
                    prev_fraction = false;
                    
                    % Store centroid to check for exponents
                    prev_centroid_y_coord = chars(i).centroid(2);
                    
                    if i+1<= num_chars
                        dist_to_next = boxes(1,i+1)-ur_x_coord;
                    else
                        dist_to_next = NaN;
                    end
                    detected = help_get_char(detected,i,chars,dist_to_next);
                end
            else
                % Normal op
                prev_fraction = false;
                
                % Store centroid to check for exponents
                prev_centroid_y_coord = chars(i).centroid(2);
                
                if i+1<= num_chars
                    dist_to_next = boxes(1,i+1)-ur_x_coord;
                else
                    dist_to_next = NaN;
                end
                detected = help_get_char(detected,i,chars,dist_to_next);
            end
            eq_string = [eq_string detected];
        case 1
            % Make sure the overlap is the next char
%             assert(overlap_idx(i+1))
            % Check if overlap is trivial
            overlap_ul = boxes(1,i+1);
            if abs(overlap_ul-ur_x_coord) <= 1
                % Normal op
                if any(strcmp(detected,control))
                    detected = ['\' detected ' '];
                end
                
                % Insert space if appropriate
                if i+1 <= num_chars && help_is_letter(detected) && help_is_letter(chars(i+1).char)
                    dist_to_next = boxes(1,i+1)-ur_x_coord;
                    if dist_to_next >= space
                        detected = [detected '\,'];
                    end
                end
                eq_string = [eq_string detected];
                
                % Store centroid to check for exponents
                prev_centroid_y_coord = chars(i).centroid(2);
            else
                overlap_char = chars(i+1).char;
                % Manually check for cases
                % If just a two '-', it is an equals
                if strcmp(detected,'-') && strcmp(overlap_char,'-')
                    eq_string = [eq_string '='];
                    
                    % Store average y_coord between the 2 equal bars
                    prev_centroid_y_coord = 1/2*(chars(i).centroid(2)+chars(i+1).centroid(2));
                end
                i = i+1; % Skip next char
            end
        otherwise
            % Grab the overlapped characters and boxes, look for bars and
            % limit-controls
            overlap_idx(i) = true;
            overlap_chars = chars(overlap_idx);
            overlap_centroids = zeros(2,length(overlap_chars));
            overlap_boxes = zeros(4,length(overlap_chars));
            bar_idx = [];
            bar_width = [];
            limit_idx = [];
            limit_height = [];
            for overlap_i = 1:length(overlap_chars)
                overlap_centroids(:,overlap_i) = overlap_chars(overlap_i).centroid;
                overlap_boxes(:,overlap_i) = overlap_chars(overlap_i).boundingbox;
                if strcmp(overlap_chars(overlap_i).char,'-')
                    bar_idx(end+1) = overlap_i;
                    bar_width(end+1) = overlap_chars(overlap_i).boundingbox(3);
                elseif any(strcmp(overlap_chars(overlap_i).char,limit_control))
                    limit_idx(end+1) = overlap_i;
                    limit_height(end+1) = overlap_chars(overlap_i).boundingbox(4);
                end
            end
            
            ul_list = overlap_boxes(1,:);
            ur_list = overlap_boxes(1,:)+overlap_boxes(3,:);
            total_width = max(ur_list)-min(ul_list);
            
            % If there is a bar that is about as wide as the entire overlap
            % region, then it is a dominant fraction bar.
            bar_idx = bar_idx(bar_width >= 0.8*total_width);
            
            % If we have a single fraction bar recognized
            if length(bar_idx) == 1
                prev_fraction = true;
                % get the bar
                frac_bar = overlap_chars(bar_idx);
                frac_y_coord = frac_bar.centroid(2);
                frac_bar_height = frac_bar.boundingbox(4);
                
                % Check if fraction is an exponent (need to check if is NaN
                % in case fraction is first character)
                if ~isnan(prev_centroid_y_coord) && frac_y_coord < prev_centroid_y_coord - 7*frac_bar_height
                    is_super = true;
                    eq_string = [strtrim(eq_string) '^{'];
                end
                
                % Store frac_bar centroid y for super/subscript checks
                prev_centroid_y_coord = frac_y_coord;
                
                % Get the numerator and denominator characters by comparing
                % centroid to the frac bar y coord.
                numer_idx = overlap_centroids < frac_y_coord;
                numer_idx = numer_idx(2,:);
                denom_idx = overlap_centroids > frac_y_coord;
                denom_idx = denom_idx(2,:);
                
                numer_eq_struct.filename = '';
                numer_eq_struct.characters = overlap_chars(numer_idx);
                denom_eq_struct.filename = '';
                denom_eq_struct.characters = overlap_chars(denom_idx);
                
                % Get subequation string
                numer_str = fn_assemble_eq(numer_eq_struct);
                denom_str = fn_assemble_eq(denom_eq_struct);
                
                detected = ['\frac{' numer_str '}{' denom_str '}'];
                
                % Limitation: fraction must be end of exponent
                if is_super
                    is_super = false;
                    detected = [detected '}'];
                end
                
                % Set counter to next character that is not in fraction
                i = find(overlap_idx,1,'last');
                
            elseif ~isempty(limit_idx)
                [~,dom_idx] = max(limit_height); % Choose largest to be the limit_char
                dom_idx = limit_idx(dom_idx);
                % Get the char
                limit_char = overlap_chars(dom_idx);
                limit_y_coord = limit_char.centroid(2);
                
                % Store limit character centroid y for super/subscript
                prev_centroid_y_coord = limit_y_coord;
                
                % Get the top and bottom characters by comparing
                % centroid to the frac bar y coord.
                top_idx = overlap_centroids < limit_y_coord;
                top_idx = top_idx(2,:);
                bottom_idx = overlap_centroids > limit_y_coord;
                bottom_idx = bottom_idx(2,:);
                
                detected = ['\' limit_char.char '\limits'];
                if any(bottom_idx)
                    bottom_eq_struct.filename = '';
                    bottom_eq_struct.characters = overlap_chars(bottom_idx);
                    bottom_str = fn_assemble_eq(bottom_eq_struct);
                    detected = [detected '_{' bottom_str '}'];
                end
                if any(top_idx)
                    top_eq_struct.filename = '';
                    top_eq_struct.characters = overlap_chars(top_idx);
                    top_str = fn_assemble_eq(top_eq_struct);
                    detected = [detected '^{' top_str '}'];
                end
                
                % Set counter to next character
                i = find(overlap_idx,1,'last');
            end
            
            eq_string = [eq_string detected];
    end
    i = i+1;
end

%% Clean up and helper functions
% If we ended in a superscript or subscript, we need to end.
if is_super || is_sub
    eq_string = [strtrim(eq_string) '}'];
end

    function tex_char = help_get_char(detected,index, chars,dist_to_next)
    % Helps recognize control characters
        if any(strcmp(detected,control)) || any(strcmp(detected,limit_control))
            detected = ['\' detected ' '];
        end
        
        % Insert space if needed between letters
        if ~isnan(dist_to_next) && help_is_letter(detected) && help_is_letter(chars(index+1).char)
            if dist_to_next >= space
                detected = [detected '\,'];
            end
        end
        tex_char = detected;
    end


end

function control = help_create_control()
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

function new_chars = help_preprocess_chars(chars,boxes)
% Special case: lim
i = 1;
new_chars(1) = chars(1);
while i <= length(chars)
    detected = chars(i).char;
    
    % If 'l' Check if 'l', '.', 'm' are close by afterwards
    if strcmp(detected,'l')
        start_x_coord = chars(i).boundingbox(1);
        % the 'm' in 'lim' should begin at around (2.333)(width of 'l')
        % after the 'l'
        end_x_coord = chars(i).boundingbox(3)*3 + start_x_coord;
        following_idx = find(boxes(1,:)>start_x_coord & boxes(1,:) < end_x_coord);
        following = chars(following_idx);
        lim_chars_idx =  (strcmp({following.char},'l')) | ...
            (strcmp({following.char},'.')) | ...
            (strcmp({following.char},'m'));
        % Translate back into chars indexing
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

