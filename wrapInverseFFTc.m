function [fh] = wrapInverseFFTc(b) 

lm = size(b,2)-1;

a = complex(zeros(2*lm,1));

% for p=0:lm 
    a(1:lm+1)=complex(b);
% end 
a = ifft(a);
%for p=1:lm+1 
fh = 2*lm*real(a(1:(lm+1)));
%end 

end