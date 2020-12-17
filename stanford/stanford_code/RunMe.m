% RunMe.m
% Jim Brewer and James Sun
% EE368 - Stanford University, Autumn 2015-16
% Final Project - LaTeX Generation from Equation Photograph
% Dec 4, 2015

%% Run this file to execute the LaTeX Generation Code on a specified image.
% Follow the comments below for directions to run the code.

clear all; close all;

%% Specify the name of the file to run the code on. 
%   The function fn_main assumes the files or folders will be in a 
%   subfolder called "Equations." To change this, modify fn_main.m

% For a clean PDF export image (comment out other):
file = 'Clean/eq1_h.jpg';

% For a real-world equation photograph (comment out other):
% file = 'Images/eq1_hr.jpg';

%% Runs the entire pipeline on the input image. Outputs a *.tex file in the
% input image directory with the same name as the input image.
%   Remove the "true" or set to false to process without showing intermediate
%   figures.
%   Add an optional third argument string to specify the name of the output
%   file (without *.tex extension):
%       "test" will place test.tex in the "Equations" folder.
%       "../test" will place test.tex in the "Code_Brewer_Sun" folder
fn_main(file, true, '../test');

