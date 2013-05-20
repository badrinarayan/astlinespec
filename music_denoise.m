function out =  music_denoise(observed, k)
%MUSIC_DENOISE Given a line spectral signal and number of frequencies
% MUSIC_DENOISE uses the rootmusic algorithm to estimate the frequencies
% and returns the debiased amplitudes and the denoised signal.
%
% out = MUSIC_DENOISE(observed, num_freqs)
% out is a struct with the following properties
% - debiased = debiased signal
% - freqs    = estimated frequencies
% - amps     = estimated amplitudes
%
% See also ROOTMUSIC, AST_DENOISE, CADZOW_DENOISE LINESPECTRUM

n = length(observed);
fs_music = sort(mod(rootmusic(observed,k),2*pi));
music_basis = exp( 1i*(0:(n-1))'*reshape(fs_music,1,k) );
fs = fs_music/2/pi;
cs = music_basis\observed;
x_music = music_basis*cs;
out = struct('debiased',x_music, 'freqs', fs, 'amps', cs);
end
