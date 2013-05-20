function out = ast_denoise(observed,varargin)
%AST_DENOISE   Denoise line spectral signal using AST and return
%frequencies and amplitudes
%
% Examples:
%
% * Simplest invocation uses only the observation
%
%   out = ast_denoise(observed)
%
% * You could pass in an estimate of the noise standard deviation to
% aid the computation of the tradeoff parameter for AST.
%
%   out = ast_denoise(observed,'noise_std',noise_std)
%
% * Alternatively, you can specify the tradeoff parameter directly
%
%   out = ast_denoise(observed,'tau',tau)
%
% By default, the backend algorithm is ADMM. You could pass in other
% options for backend like SDPT3 and CVX.
%
%   out = ast_denoise(observed,'backend','admm')  % (default)
%   out = ast_denoise(observed,'backend','sdpt3') % (SDPT3 backend)
%   out = ast_denoise(observed,'backend','cvx')   % (CVX backend)
%
% The Lasso backend produces an approximate solution and can be used 
% in conjunction with a grid size parameter (default 2^12)
% 
% Given samples from a mixture of exponentials and an estimate of the
% standard deviation of noise, this code runs the AST algorithm to find
% the frequencies and amplitudes of the mixture of exponentials
%   
% It returns a structure out with the following outputs
%
% Tradeoff Parameter:
%   tau          - Tradeoff parameter used
% Denoised Signal:
%   out.estimate - estimated signal
%   out.debiased - debiased signal using the estimated frequencies
% Estimated Frequencies and Amplitudes:
%   out.freqs    - estimated frequencies (using dual polynomial)
%   out.amps     - estimated amplitudes (debiased)
%
% Dual Polynomial (only for AST, CVX and SDPT3)
%
% out.pts      - frequency values for dual polynomial
% out.polyval  - dual polynomial values
%
% See also AST_ADMM, AST_SDPT3, AST_CVX, MUSIC_DENOISE, CADZOW_DENOISE

p = inputParser;
addRequired(p,'observed');
n = length(observed);
addParamValue(p,'backend','admm',...
    @(x) any(validatestring(x,{'admm','cvx','sdpt3','lasso'})));
addParamValue(p,'noise_std',1.5*estimate_noise_std(observed));
addParamValue(p,'gridsize',2^12);
addParamValue(p,'tau',0);
parse(p,observed,varargin{:});

noise_std = p.Results.noise_std;
tau = p.Results.tau;
if(tau==0)
    tau = sqrt(log(n)+log(4*pi*log(n)))*noise_std;
end

pts = [];
polyval = [];
Tu = [];

switch p.Results.backend
    case 'admm'
        [x_ast,dual_poly_coeffs,Tu] = ast_admm(observed,tau);
        [debiased,fs,cs,pts,polyval] = dual_poly_debias(dual_poly_coeffs,observed);
    case 'sdpt3'
        [x_ast,dual_poly_coeffs,Tu] = ast_sdpt3(observed,tau);
        [debiased,fs,cs,pts,polyval] = dual_poly_debias(dual_poly_coeffs,observed);
    case 'cvx'
        [x_ast,dual_poly_coeffs,Tu] = ast_cvx(observed,tau);
        [debiased,fs,cs,pts,polyval] = dual_poly_debias(dual_poly_coeffs,observed);
    case 'lasso'
        [debiased,fs,cs,x_ast] = ast_lasso(observed,sqrt(n)*tau,p.Results.gridsize);
    otherwise
        disp('Unknown backend. The following work: admm,sdpt3,cvx');
        exit
end

out = struct('freqs',fs(:),'amps',cs(:),'debiased',debiased,...
    'estimate',x_ast,'Tu',Tu,'pts',pts,'polyval',polyval,...
    'tau',tau,'noise_std',noise_std);
end
