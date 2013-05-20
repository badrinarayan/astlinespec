function [xhat,dual_poly_coeffs,Tu] =  ast_cvx(y,tau)
% AST_CVX Solves AST for line spectral signals using CVX
%
% [xhat,dual_poly_coeffs,Tu] =  AST_SDPT3(y,tau)
%
% Denoises uniform samples of a mixture of complex sinusoids and
% returns the solution x of the atomic soft thresholding problem
%            1            2
% minimize  ---|| y - x || + tau ||x|| (AST)
%    x       2                        A
% where the atoms are complex sinusoids. AST is recast as an SDP
% solved using SDPT3.
%
% The function returns the optimal value x, and the dual polynomial coefficients
% q = (y - x)/mu. The dual polynomial coefficients specify a trigonometric
% polynomial whose absolute value reaches unity at the supporting
% frequencies of the solution.
% 
% This function requires the installation of CVX package which is
% available at http://cvxr.com/cvx/
%
% See also AST_SDPT3, AST_ADMM

n = length(y);
cvx_begin quiet
variable u0;
variable u(n-1) complex;
variable Tu(n,n) hermitian;
variable t;
variable xhat(n) complex;
minimize (tau*(n*u0 + t)/2 + 1/2*pow_pos(norm(y - xhat,'fro'),2))
subject to
Tu == toeplitz([u0;conj(u)],[u0;conj(u)]');
[Tu xhat;
 xhat' t] == semidefinite(n+1,n+1);
cvx_end
dual_poly_coeffs = (y - xhat)/tau;
end