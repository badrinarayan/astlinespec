function out = m1func(~,f0s,c1s,f1s,n)
%M1FUNC Compute metric m1 = Sum of absolute amplitudes in far region
%
%   Inputs:
%
%   c0s = original amplitudes
%   f0s = original frequencies
%   c1s = estimated amplitudes
%   f1s = estimated frequencies
%   n = number of samples
%
%   Outputs:
%
%   m1 Sum of absolute amplitudes in far region
%   See Theorem 2 of http://arxiv.org/abs/1303.4348v1
k = length(f0s);

near_amps = 0;
for j=1:k
    f0 = f0s(j);    
    near_amps = near_amps + sum(abs(c1s(dtorus(f1s,f0)<0.16/n)));
end
out = sum(abs(c1s)) - near_amps;
end