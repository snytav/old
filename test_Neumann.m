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

u = neumannsolver(hx,hz,left,cent,righ,fr1,fr2,fz1,fz2,b);

uf = dlmread('uf.txt');
u1 = reshape(uf,im+2,lm+1);
u2 = reshape(u',(im+2)*(lm+1),1);

[m,i] = max(abs(u2-uf));

qq = 0;


