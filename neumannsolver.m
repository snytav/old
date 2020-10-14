%     hx,hz          - radial and "z" grid steps
%     left,cent,righ - radial stencil
%     fr1,fr2        - radial (Dirichlet) boundary conditions
%     fz1,fz2        - "z" (Neumann) boundary conditions
%     b              - right-hand side
%     u              - solution 
function [u] = neumannsolver()%(hx,hz,left,cent,righ,fr1,fr2,fz1,fz2,b) 


b = dlmread('b_3.txt');

fr1 = dlmread('fr1.txt');
fr2 = dlmread('fr2.txt');

fz1 = dlmread('fz1.txt');
fz2 = dlmread('fz2.txt');
lm = size(fr1,1)-1;
im = size(fz1,1)-2;

b = reshape(b,im+2,lm+1);

a = dlmread('hxhz.txt');
hx = a(2);
hz = a(3);

left = dlmread('l.txt');
cent = dlmread('c.txt');
righ = dlmread('r.txt');
hz2 = hz*hz;
hx2 = hx^2;

for i = 0:lm
    lambda(i+1) = 4/hz2*sin(i*pi/2/lm)^2; % !i!dcos(pi/lm*i)!1d0
end

for i = 1:im+2
    b(i,1)    = b(i,1)    - 2*fz1(i)*hx2;
    b(i,lm+1) = b(i,lm+1) + 2*fz2(i)*hx2; 
end

for k = 0:lm
	b(2,k+1)    = b(2,k+1)    - left(2)*fr1(k+1)*hx2;
	b(im+1,k+1) = b(im+1,k+1) - righ(im+1)*fr2(k+1)*hx2;
end

left1 = zeros(im+2,1);
righ1 = zeros(im+2,1);
fh = zeros(im+2,lm+1);
for i = 2:im+1
	fh(i,:) = directFT(b,i);
	left1(i) = left(i)*hx2;
	righ1(i) = righ(i)*hx2;
end

for k = 0:lm
	for i = 1:im+2
   	    cent1(i) = hx2*(cent(i) + lambda(k+1));
    end

    f = fh(:,k+1);
    x = crprogon(left1,cent1,righ1,f);
    u(:,k+1) = x;
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end 

f = zeros(im+2,lm+1);

for i = 2:im+1
	fhw = inverseFT(u,i);
    f(i,:) = fhw;
end
u = f*(2/lm);
    
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
for k = 1,lm+1
	u(1,k)    = fr1(k);
	u(im+2,k) = fr2(k);
end

end 