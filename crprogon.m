function [x] = crprogon(a,b,c,f) 

im = size(f,1)-2;

be(1)=f(1)/b(2);
al(1)=-c(1)/b(2);
al = zeros(im+2,1);
be = zeros(im+2,1);
for i=2:im+2
   if i  == im+2
       qq = 0;
   end
   al(i)  = -c(i)/(b(i-1)+al(i-1)*a(i));
   be(i)  = (f(i)-a(i)*be(i-1))/(b(i-1)+al(i-1)*a(i));
end


x(im+2) = be(im+2);
for i=im+1:-1:2 
    x(i)=al(i)*x(i+1)+be(i);
end 

end