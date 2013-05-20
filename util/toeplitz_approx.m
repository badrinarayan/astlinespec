function Ta = toeplitz_approx(A)
% TOEPLITZ_APPROX(A)
% For a square matrix A, find the Toeplitz approximation
% by averaging along the diagonals.
[M,N] = size(A);
Ta = zeros(size(A));
b = zeros(M+N-1,1);
for n = -(M-1):(N-1)
		b(M+n) = mean(diag(A,n));
end
Ta = T(b);
end
