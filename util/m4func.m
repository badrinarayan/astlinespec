function out = m4func(~,f0s,~,f1s,n)
%M3FUNC Compute metric m3 = Max Deviation of frequencies
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
%   m3 Max difference between true amplitude and sum of estimated amplitudes
%   in the near region
%   See Theorem 2 of http://arxiv.org/abs/1303.4348v1
k = length(f0s);

out = 0;
for j=1:k
    out = max(out,min(dtorus(f0s(j),f1s)));
end

end