Jim Brewer and James Sun
EE368 - Stanford University, Autumn 2015-16
Final Project - LaTeX Generation from Equation Photograph
Dec 4, 2015

LaTeX Generation:
To run the basic pipeline process, run the RunMe.m MatLab script.
Inside this script you can specify which file you would like to process.
See script for further details.


Test Files:
Included in the Equations/ folder are sample images of equations the
code generation pipeline can be run on. This includes "clean" images 
extracted directly from PDF files of equations and real photographs
of equations with uneven lighting and skew.


Demo:
In the Demo/ folder is a batch file "mat.bat" that can be run to grab any images in
the Demo/Images directory and batch process them. It will export the *.tex file and
compile the *.tex file and show the final PDF if a LaTeX is installed on the computer.


Test Scripts:
Included are scripts that were used for testing and generating test data for the report.

AssembleTest.m
CombinedOptTest.m
fn_lighting_otsu.m
fn_resizeForTest.m
fn_rotateForTest.m
LightingTest.m
LightingTimingTest.m
SkewTest.m
testEquationMatching.m
testMatching.m

Test Data Files:
equation_truth.mat
equation_truth_mod.mat
EquationTestData.mat