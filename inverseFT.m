function [fhw] = inverseFT(b,i)

lm = size(b,2) - 1;

for k = 0:lm
	fh(i,k+1) = 0

	if k == 0 || k == lm
	   rho(k+1) = 0.5*b(i,k+1)
    else
	   rho(k+1) = b(i,k+1)
    end
end

fhw = wrapInverseFFTc(rho);


end