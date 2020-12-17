% Tests fn_lighting_compensation. Expects equation photographs to be in the
% "Lighting" directory.
addpath('../');
lighting_dir = 'Lighting\';
results_dir = [lighting_dir 'Results\'];
listing = dir(lighting_dir);

for i = 1:length(listing)
    fname = listing(i).name;
    if ~listing(i).isdir && ~strcmp(fname,'.') && ~strcmp(fname,'..')
        img = imread([lighting_dir fname]);
        bw_img = fn_lighting_compensation(img);
        imwrite(bw_img,[results_dir 'result_' fname]);
    end
end

