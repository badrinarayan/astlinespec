%% AST Experiments
% 
% Here we show how to replicate the experimental setup in
%
% [1] Badri Narayan Bhaskar, Gongguo Tang, Benjamin Recht, 
% "Atomic norm denoising with applications to line spectral estimation"

%% Parameters
% The experiments in one runs 20 iterations for these combination
% of parameters (1260 in total). We demonstrate how to run a
% random experiment from this collection

ns = [64,128,256];     % number of samples
k_factors = [4,8,16];  % n/k, k = number of frequencies
SNRs = -10:5:20;       % signal to noise ratio

%% Pick a random experiment
pickone  = @(xs) xs(ceil(length(xs)*rand));

n        = pickone(ns);
k_factor = pickone(k_factors);
k        = n/k_factor;
SNR      = pickone(SNRs);

% Or make your own
n = 64;
k = 8;
SNR = 15;

fprintf('n = %d, k = %d, SNR = %d\n',n,k,SNR)

%% Generate the line spectral signal with these parameters

[signal,amps,freqs] = linespectrum(n,k,...
            'amplitudes',randn(k,1).^2.*exp(1i*2*pi*rand(k,1)),...
            'frequency_spacing','random',...
            'minimum_spacing',1/n);
noise_std = norm(signal)/sqrt(n)*10^(-SNR/20);
observed = signal + noise_std*(randn(n,1) + 1i*randn(n,1))/sqrt(2);

%% Run AST, MUSIC, Cadzow and Lasso to denoise

fprintf('Running AST...')
tic;ast     = ast_denoise(observed);
fprintf('Finished in %f s.\n', toc);

fprintf('Running Cadzow...')
tic;cadzow  = cadzow_denoise(observed,k);
fprintf('Finished in %f s\n', toc);

fprintf('Running MUSIC...')
tic;music   = music_denoise(observed,k);
fprintf('MUSIC finished in %f s.\n', toc);

fprintf('Running Lasso...')
tic;lasso12   = ast_denoise(observed,'backend','lasso');
fprintf('Lasso finished in %f s.\n', toc);


%% Compute Mean Square Errors and Frequency Deviation Metrics

mse_func = @(estimate) norm(estimate.debiased(:)-signal(:))/norm(signal);
m1       = @(estimate) m1func(estimate.amps,estimate.freqs,amps,freqs,n);
m2       = @(estimate) m2func(estimate.amps,estimate.freqs,amps,freqs,n);
m3       = @(estimate) m3func(estimate.amps,estimate.freqs,amps,freqs,n);

fprintf('\nAST    MSE = %.4f\n', mse_func(ast));
fprintf('MUSIC  MSE = %.4f\n', mse_func(music));
fprintf('Cadzow MSE = %.4f\n', mse_func(cadzow));
fprintf('Lasso MSE = %.4f\n\n', mse_func(lasso12));

% Frequency Metrics
fprintf('Frequency Metrics\n');
fprintf('\nm1:\n AST    = %g\n MUSIC  = %g\n Cadzow = %g\n Lasso  = %g',...
    m1(ast),m1(music),m1(cadzow),m1(lasso12));
fprintf('\nm2:\n AST    = %g\n MUSIC  = %g\n Cadzow = %g\n Lasso  = %g',...
    m2(ast),m2(music),m2(cadzow),m2(lasso12));
fprintf('\nm3:\n AST    = %g\n MUSIC  = %g\n Cadzow = %g\n Lasso  = %g\n',...
    m3(ast),m3(music),m3(cadzow),m3(lasso12));