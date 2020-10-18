function [f] = inverseFT2D(u)
 [im2,lm1] = size(u);
 im = im2-2;
 lm = lm1-1;
 f = zeros(im+2,lm+1);
 
 for i = 2:im+1
 	fhw = inverseFT(u,i);
    f(i,:) = fhw;
 end
end
