function [ name ] = writeTex(eq_string, fileName )
%writeTex Write the passed equation string into a compilable .tex file
%   with passed filename.

fileID = fopen(strcat(fileName,'.tex'),'w');

header = {'\\documentclass[20pt]{report}\n'
    '\\DeclareMathSizes{12}{30}{30}{30}\n'
    '\\begin{document}\n'
    '\\begin{center}'
    '\n\\(\n'};

for i = 1:length(header)
    fprintf(fileID, header{i});
end

% Replace \ with with \\
fprintf(fileID,strrep(eq_string,'\','\\'));

footer = {'\n\\)\n'
    '\\end{center}\n'
    '\\end{document}\n'};
for i = 1:length(footer)
    fprintf(fileID, footer{i});
end

fclose(fileID);

end

