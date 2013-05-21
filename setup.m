%SETUP Installation
javaaddpath java
addpath redir
addpath util
addpath backends
builddocsearchdb('html')
disp('Congratulations! You have now successfully installed Atomic Line Spectral Estimation Toolbox. ');
disp('cd to redir directory and install CVX and SDPT3 if you need those backends.');
disp('To get started, type: showdemo astdemo');
disp('For a general overview, type help astlinespec');

