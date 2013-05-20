function [ poles, amplitudes ] = poles_amps( signal, num_poles )
%POLES_AMPS Get the poles and amplitudes from the signal
%   Using the Matrix Pencil algorithm, determine the
% location of the poles and the amplitudes from the signal
signal = signal(:);
n = length(signal);
z = prony_pencil_est(signal,num_poles);
w = mod(2*pi+atan2(imag(z),real(z)),2*pi).';
c = exp(1i*w*(0:n-1)).'\signal;
[~,I] = sort(w);
poles = w(I)/2/pi;
amplitudes = c(I);
end

