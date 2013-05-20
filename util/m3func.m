function out = m3func(c0s,f0s,c1s,f1s,n)
%M3FUNC Compute metric m3 = Max Deviation of amplitudes in near region
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
    out = max(out, abs(c0s(j) - sum(c1s(dtorus(f1s,f0s(j))<0.16/n))));
end

end