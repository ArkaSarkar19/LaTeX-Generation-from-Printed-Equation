% Tests Equation assembly on the 10 test equation segmentations found in
% EquationTestData.mat
addpath('../');
load EquationTestData
for i = 1:length(equations)
    A = fn_assemble_eq(equations(i));
    disp([num2str(i) ' => ' A])
end