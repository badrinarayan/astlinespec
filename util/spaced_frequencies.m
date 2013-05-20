function [fs,success,iteration] = spaced_frequencies(separation, m)
%SPACED_FREQUENCIES returns well sepearated frequencies in [0,1]
%
% [fs,success,iteration] = spaced_frequencies(m, k)
% 
% Given a minimum separation m in [0,1] and a number of frequencies k,
% SPACED_FREQUENCIES returns k frequencies that are separated by atleast m.
%
% The algorithm procedes iteratively and backs off and rejects when the
% condition is violated. It optionally returns a binary success flag and 
% the number of iterations needed to return a valid random collection.
%
% See also LINESPECTRUM

fs = zeros(m,1);
max_iterations = 10000;
iteration = 1;
while iteration < max_iterations
  r = RandSep(separation);
  success = 1;
  for i=1:m
    fs(i) = r.next();
    if fs(i)==-1
      success = 0;
      break;
    end
  end
  if success
    break;
  end
  iteration = iteration + 1;
end
end
