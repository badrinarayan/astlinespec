function z = prony_pencil_est(x,k)
% PRONY_PENCIL_EST Find the complex poles by Matrix Pencil algorithm
%
% For details see
%
% Yingbo Hua, Tapan Sarkar, "Matrix Pencil Method for Estimating Parameters
% of exponentially Damped/Undamped Sinusoids in Noise", IEEE Transactions on
% Acoustics, Speech and Signal Processing, Vol. 38, No. 5, May 1990
%
n = length(x);
X0 = hankel(x(1:n-k),x(n-k:n-1));
X1 = hankel(x(2:n-k+1),x(n-k+1:n));
[V,D] = eig(pinv(X0)*X1);
z = diag(D).';
end