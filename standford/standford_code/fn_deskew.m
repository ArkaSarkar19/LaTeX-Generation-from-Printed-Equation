function [ deskew_img ] = fn_deskew( bw_img )
%FN_DESKEW NOT FOR USE. USE fn_deskew2 instead.Attempts to implement a fast
%skew detection algorithm using the frequency space. Results are promising 
%but not good enough to supplant Hough transform method.


F=fft2(bw_img);
[H,W] = size(F);

S=fftshift(F);
A=log(1+abs(S));
% Simulate a high pass filter
center_perc = 0.1;
A(floor(H/2-center_perc*H):floor(H/2+center_perc*H),...
    floor(W/2-center_perc*W):floor(W/2+center_perc*W)) = 0;

% Divide into quadrants
quad_h = floor(H/2)+mod(H,2);
quad_w = floor(W/2)+mod(W,2);

quadrants = zeros(quad_h,quad_w, 4);
for i = 1:4
    m = mod(i+1,2)+1;
    n = floor((i+1)/2);
    % Duplicate center axes
    if m > 1
        m_int = (m-1)*quad_h:m*quad_h-1;
    else
        m_int = (m-1)*quad_h+1:m*quad_h;        
    end
    if n > 1
        n_int = (n-1)*quad_w:n*quad_w-1;
    else
        n_int = (n-1)*quad_w+1:n*quad_w;        
    end
    quad = A(m_int,n_int);
    quadrants(:,:,i) = quad;
end

angles = zeros(1,2);
for i = 1:4
    quad = quadrants(:,:,i);
    sorted = sort(quad(:),'descend');
    bw_quad = quad > sorted(20);
    [y,x] = find(bw_quad);
    [h,w] = size(bw_quad);
    angles(i) = help_ls(x,y,h,w,i);
%     switch i
%         case 1            
%             bw_quad(end,end) = 1;
%         case 2
%             bw_quad(1,end) = 1;
%         case 3
%             bw_quad(end,1) = 1;
%         case 4
%             bw_quad(1,1) = 1;
%     end
%             
%     % Limit theta to -45 to 45;
%     [H,T,~] = hough(bw_quad,'RhoResolution',1,'Theta', -45:0.1:45);
%     P = houghpeaks(H);
%     if i == 1 || i == 4
%         angles(i) = T(P(2));    
%     else
%         angles(i) = -T(P(2));
%     end
end
% detected_angle = mean(angles);
angles(isnan(angles)) = [];
if isempty(angles)
   detected_angle = 0; 
else
    var_limit = var([0 45]);
    if var(angles) >= var_limit
%         detected_angle = angles(1);
        detected_angle = mean(angles);
    else
        detected_angle = mean(angles);
    end
    if detected_angle >= 30
        detected_angle = 0;
    end
end

deskew_img = detected_angle;

end

function theta = help_ls(x,y,h,w,quad)
    % Compensate for 1-indexing
    x = x-1;
    y = y-1;
    y = -y;
    
    % x*cos + y*sin = rho
    % [x y][cos;sin] = rho
    % [x y][cos/rho; sin/rho] = 1;
    % (sin/rho)/(cos/rho) = sin/cos = tan
    A = [x y];
    origin_weight = 15;
    switch quad
        case 1
            added = repmat([w,-h],origin_weight,1);
        case 2
            added = repmat([w,0],origin_weight,1);
        case 3
            added = repmat([0,-h],origin_weight,1);
        case 4
            added = repmat([0,0],origin_weight,1);
    end
    A = [A;added];
    
    b = ones(size(A,1),1);
    trigs = (A'*A)^-1*A'*b;
    theta = atan2d(trigs(2),trigs(1));
    
%     A = [x -y ones(size(y))];
%     [~,~,v] = svd(A);
%     v_small = v(:,end);
%     theta = atand(v_small(2)/v_small(1));
%     rho = v_small(3);
%     y_prime = (rho-v_small(1)*x)/v_small(2);
    J = sum((A*trigs-b).^2);
    if J >= 0.1
        theta = NaN;
    else
        switch quad
            case 1
                theta = theta;
            case 2
                theta = 90+theta;
            case 3
                theta = 90+theta;
            case 4
                theta = -(90+theta);
        end
    end
 
        
end

