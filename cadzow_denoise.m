function out = cadzow_denoise(y, r)
%CADZOW_DENOISE Given a line spectral signal and number of frequencies
% CADZOW_DENOISE uses Cadzow's algorithm to estimate the frequencies
% and returns the debiased amplitudes and the denoised signal.
%
% out = CADZOW_DENOISE(observed, num_freqs)
% out is a struct with the following properties
% - debiased = debiased signal
% - freqs    = estimated frequencies
% - amps     = estimated amplitudes
%
% Cadzow's algorithm alternately projects the associated toeplitz matrix of
% the noisy line spectral signal onto Rank-num_freqs Grasmmanian manifold
% and space of Toeplitz matrices. The frequencies are then extracted using
% matrix pencil algorithm. The debiased amplitudes and the reconstructed 
% signal are returned.
%
% References
%
% James Cadzow, "Signal enhancement - a composite property mapping 
% algorithm," IEEE Trans. on Acoustics, Speech and Signal Processing,
% vol. 36, no. 1, pp. 49?62, 1988.
%
% Yingbo Hua, Tapan Sarkar, "Matrix Pencil Method for Estimating Parameters
% of exponentially Damped/Undamped Sinusoids in Noise", IEEE Transactions on
% Acoustics, Speech and Signal Processing, Vol. 38, No. 5, May 1990
%
% See also AST_DENOISE, MUSIC_DENOISE, LINESPECTRUM, MATRIX_PENCIL_DEBIAS

x = y;
max_iterations = 5000;
tol = 1e-10;
ratio = tol;

n = length(x);
X = T(x);

for iter=1:max_iterations
  if (ratio<tol), break; end
  % Alternate Projection using SVD, toeplitz_approx
  [U,D,V] = svd(X);
  V((r+1):end,(r+1):end)=0;
	d = diag(D);
  X = toeplitz_approx(U(:,1:r)*diag(d(1:r))*V(:,1:r)');
  % Update convergence parameter
  ratio = D(r+1,r+1)/D(r,r);
end

x = invT(x,n);
[x_debiased,fs,cs] = matrix_pencil_debias(invT(X,n), y, r);
out = struct('debiased',x_debiased, 'freqs', fs, 'amps', cs);
end
