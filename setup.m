%SETUP Installation
javaaddpath java
addpath redir
addpath util
addpath perf
addpath backends
builddocsearchdb('html')
disp('cd to redir directory and install CVX and SDPT3 if you need those backends.');
disp('To get started, type: showdemo astdemo');
disp('For a list of files');

