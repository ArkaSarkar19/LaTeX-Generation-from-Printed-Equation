% Tests both fn_lighting_compensation and fn_deskew2 on the photographed
% images in the "Lighting" Directory.
addpath('../');
lighting_dir = 'Lighting\';
combined_dir = 'Opt\';
results_dir = [combined_dir 'Results\'];
listing = dir(lighting_dir);

for i = 1:length(listing)
    fname = listing(i).name;
    if ~listing(i).isdir && ~strcmp(fname,'.') && ~strcmp(fname,'..')
        img = imread([lighting_dir fname]);
        bw_img = fn_lighting_compensation(img);
        [out_img, detected_angle] = fn_deskew2(bw_img,true,true,5);
        imwrite(out_img,[results_dir 'result_hough_' fname]);
    end
end
