function [fh] = wrapDirectFFTc(b) 

lm = size(b,2)-1;

a = complex(zeros(2*lm,1));

for p=0:lm 
    a(p+1)=complex(b(p+1),0);
end 
a = fft(a); 

%for p=1:lm+1 
fh = real(a(1:(lm+1)));
%end

end