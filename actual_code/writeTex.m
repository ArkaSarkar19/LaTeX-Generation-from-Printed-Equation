function [ name ] = Final_code(final_string )



header = {'\\documentclass[20pt]{report}\n'
    '\\DeclareMathSizes{12}{30}{30}{30}\n'
    '\\begin{document}\n'
    '\\begin{center}'
    '\n\\(\n'};

for i = 1:length(header)
    fprintf(fileID, header{i});
end

% Replace \ with with \\
fprintf(fileID,strrep(final_string,'\','\\'));

footer = {'\n\\)\n'
    '\\end{center}\n'
    '\\end{document}\n'};
for i = 1:length(footer)
    fprintf(fileID, footer{i});
end

fclose(fileID);

end

