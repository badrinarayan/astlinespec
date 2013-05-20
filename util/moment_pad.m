function z=moment_pad(x,nsmps,m)
z = ifft(x,nsmps)*sqrt(nsmps);
if length(z)==0 
	z = zeros(m,1);
else
	z = z(1:m);
end
