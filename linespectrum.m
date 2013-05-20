function [signal,amps,freqs] = linespectrum(num_samples,num_freqs,varargin)
%LINESPECTRUM Construct a line spectral signal
%
% [signal, amps, freqs] = LINESPECTRUM(num_samples,num_freqs)
%
% Generate n dimensional line spectral signal such that
% 
% x(m+1) = sum of amps(l)*exp(1i*2*pi*freqs(l)m)
% for m = 1 to num_samples, for l = 1 to num_freqs
%
% By default, the num_freqs frequencies are randomly chosen between 0 and 1
% with minimum spacing of 1/n. But all of these can be customized with
% optional arguments. The function also returns the chosen amplitudes and
% frequencies.
%
% signal = LINESPECTRUM(100,5,'frequency_spacing','equispaced')
%
% This invocation returns 5 equispaced frequencies between 0 and 1 and
% provides 100 uniform samples of a line spectral signal with those
% frequencies.
% 
% If the frequency spacing is random (default), you can specify a minimum
% spacing between the frequencies (which defaults to 1/num_samples) using
% the following invocation - any floating point number between 0 and 1 is a
% valid candidate.
% 
% [signal, amps, freqs] = LINESPECTRUM(num_samples,num_freqs,...
%    'minimum_spacing',0.5/num_samples)
%
% The preceding call uses SPACED_FREQUENCIES to generate freqs that respect
% the minimum separation condition. SPACED_FREQUENCIES uses rejection
% sampling to return the frequencies and may fail if it is infeasible or if
% there aren't sufficient iterations.
%
% Instead of specifying random or equispaced frequencies, it is also
% possible to manually specify them as in the following example:
%
% signal = LINESPECTRUM(100,3,'frequencies',[0.2,0.3,0.6])
%
% The amplitudes default to unity but it is possible to manually specify
% the amplitudes as in the following example which uses complex normal
% amplitudes with a magnitude of 1 (random phases).
% 
% signal = LINESPECTRUM(100,5,'amplitudes',exp(1i*2*pi*randn(5,1)))
% 
% See also AST_DENOISE, SPACED_FREQUENCIES

% Parse the input
p = inputParser;
addRequired(p,'num_samples',@isnumeric);
addRequired(p,'num_freqs',@(k) isnumeric(k) && (k<num_samples/2));
addParamValue(p,'frequency_spacing','random',...
    @(x) any(validatestring(x,{'equispaced','random'})));
addParamValue(p,'minimum_spacing',1/num_samples,@isnumeric);
addParamValue(p,'amplitudes',ones(num_freqs,1));
addParamValue(p,'frequencies',-ones(num_freqs,1),@(x) (length(x)==num_freqs));
parse(p,num_samples,num_freqs,varargin{:});

amps = p.Results.amplitudes;
freqs = p.Results.frequencies;

if sum(freqs)+num_freqs==0 % Means Default
  if strcmp(p.Results.frequency_spacing,'equispaced')
    freqs = mod(linspace(0,1,num_freqs+1).'+rand,1);
    freqs(end) = [];
  else
    freqs  = spaced_frequencies(p.Results.minimum_spacing,num_freqs);
  end
end

freqs = sort(freqs);
freqs = freqs(:);
signal = sum(diag(amps)*exp(1i*2*pi*freqs*(0:num_samples-1)),1).';

end