function [fh] = directFT2D(b)
     [im2,lm1] = size(b);
     im = im2-2;
     lm= lm1-1;
     
    fh = zeros(im+2,lm+1);
    for i = 2:im+1
        fh(i,:) = directFT(b,i);
%         left1(i) = left(i)*hx2;
%         righ1(i) = righ(i)*hx2;
    end
end