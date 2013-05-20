function x = invT(Tx,n)
% Inverts the Toeplitz operator for complex entries
% x = invT(Tx,n)
% Find x so that Tx = T(x)
x = [Tx(end:-1:1,1); Tx(1,2:end).'];
end