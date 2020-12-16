SET dir=%~dp0
matlab -nosplash -nodesktop -minimize -r "run('%dir%\demo_latex_gen.m');quit"; 
