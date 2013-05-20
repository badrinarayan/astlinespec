function out = m2func(~,f0s,c1s,f1s,~)
%M2FUNC Compute m2 = Weighted mean-squared deviation of frequencies
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
%   m2 Mean-squared deviation of frequencies weighted by absolute 
%   estimated amplitudes
%   See Theorem 2 of http://arxiv.org/abs/1303.4348v1

kh = length(c1s);

out = 0;
for j=1:kh
    out = out + abs(c1s(j))*min(dtorus(f1s(j),f0s).^2);
end

end