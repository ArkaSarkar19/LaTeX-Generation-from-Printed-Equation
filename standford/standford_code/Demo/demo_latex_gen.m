% This demo script should be in a separate directory named "Demo" that is
% one directory below the main LaTeX_Gen project Directory. Within this 
% "Demo" directory, we also assume the existance of an "Images" directory, 
% in which the equation photos are kept, and a "Finished" directory to
% which processed photos can be moved. 
% Also, this script was run on a Windows 10 machine with the following 
% components installed:
%   1. a LaTeX compiler
%   2. a PDF Viewer
%   3. Matlab R2015b with Image Processing Toolbox
%   4. Wordpad
% The authors cannot guarantee proper function of this script with all 
%   system configurations.


clear;
close all;

my_full_path = which('demo_latex_gen');
[file_dir,~,~] = fileparts(my_full_path);
cd(file_dir);

box_photo_dir = fullfile(file_dir,'Images\');

addpath('..');
listing = dir(box_photo_dir);

for i = 1:length(listing)
    fname = listing(i).name;
    [~,~,ext] = fileparts(fname);
    if strcmp(ext,'.jpg')
        [eq_string, output_file] = fn_demo(fullfile(box_photo_dir,fname),true);
        disp(eq_string);
        system(['pdflatex ' output_file]);
        winopen([output_file '.pdf']);
        system(['write ' output_file '.tex']);
        
        movefile(fullfile(box_photo_dir,fname) ,fullfile(file_dir,'Finished',fname));
        delete('*.aux');
        delete('*.log');
		if i+1 <= length(listing)
			h = helpdlg(['Continue to ' listing(i+1).name '?'], 'Continue?');
			waitfor(h);
		end

    end
end


h = helpdlg('Ready to exit?', 'Exit?');
waitfor(h);