% Compares the run times for our implementation of adaptive
% thresholding using mean filtering vs using otsu's method.
addpath('../');
fname = '../Equations/Images/demo_equation_hr.jpg';
img = imread(fname);

t = cputime;
tic
for i = 1:20
    fn_lighting_compensation(img);
end
toc
e1 = cputime-t;

t = cputime;
tic
for i = 1:20
    fn_lighting_otsu(img);
end
toc
e2 = cputime-t;
% These numbers seem inaccurate, but the relative times taken is still
% useful to compar.
disp(['Our implementation: ' num2str(e1)])
disp(['Using Otsu''s: '  num2str(e2)])