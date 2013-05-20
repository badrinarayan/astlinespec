Atomic Line Spectral Estimation Toolbox
========================================

---------------------------------------------------------------------------

Copyright (2013) Badri Narayan Bhaskar, Gongguo Tang, Benjamin Recht

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

---------------------------------------------------------------------------

This set of MATLAB files contain implementations of Atomic Soft 
Thresholding algorithms for denoising line spectral signals as described in
the paper 

"Atomic Norm Denoising for Line Spectral Estimation", by Badri Narayan 
Bhaskar, Gongguo Tang and Benjamin Recht (2012). Preprint available at
http://arxiv.org/abs/1204.0562

The algorithm is implemented in ast_denoise.m and scripts in the subfolder
backends. For usage details, type

help astlinespec

in the MATLAB prompt. The code also redistributes SpaRSA script for solving
the Lasso problem (See http://www.lx.it.pt/~mtf/SpaRSA/) and the CVX and 
SDPT3 packages for semidefinite programming available at http://cvxr.com and
http://www.math.nus.edu.sg/~mattohkc/sdpt3.html



