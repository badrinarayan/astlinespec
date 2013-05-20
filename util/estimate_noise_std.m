function noise_std = estimate_noise_std(observed)
%ESTIMATE_NOISE_STD Estimate the standard deviation of noise in mixture of 
% exponentials
%
% noise_std = estimate_noise_std(observed)
%
% See also LINESPECTRUM
n       = length(observed);
[~,R]   = corrmtx(observed,round(n/3),'covariance'); % R = Exx' + \sigma^2 I
all_evs = sort(eig(R),'descend');
L = length(all_evs);
noise_evs = all_evs(ceil(3*L/4):end);
% At least half the eigenvalues must correspond to noise.
noise_std = sqrt(mean(noise_evs));
end
