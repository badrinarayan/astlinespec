function [x, fs, cs] = matrix_pencil_debias( xhat, y, num_freq )
%MATRIX_PENCIL_DEBIAS Debias using Matrix Pencil Prony
%
% x = MATRIX_PENCIL_DEBIAS(xhat,y,num_freq) 
% 
% Returns the debiased signal using matrix pencil. xhat is the denoised
% signal, y is the original signal and num_freq (optional) is the number of
% frequencies

n = length(y);
debias_tol = 0.1;

Tx = T(xhat);
[p,q] = size(Tx);
p = min(p,q);
Tx = Tx(1:p,1:p); 

[~,E] = eig(Tx);
e = diag(E);
idx = (e>debias_tol*max(e)); 
if nargin==3
  [fs, ~ ] = poles_amps(xhat, num_freq);
else
  [fs, ~] = poles_amps(xhat, sum(idx) );
end
est_basis = exp( 1i*(0:(n-1))'*reshape(fs,1,num_freq) );

% fit the coefficients via least squares
cs = est_basis\y ;
% return the debiased signal
x = est_basis*cs;
end

