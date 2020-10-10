c     hx,hz          - radial and "z" grid steps
c     left,cent,righ - radial stencil
c     fr1,fr2        - radial (Dirichlet) boundary conditions
c     fz1,fz2        - "z" (Neumann) boundary conditions
c     b              - right-hand side
c     u              - solution 
	subroutine neumannsolver(hx,hz,left,cent,righ,
     =                       	 fr1,fr2,fz1,fz2,b,u) 
      implicit none
      include 'four.par'
	real*8 hz,b(im+2,0:lm),u(im+2,0:lm),pi,lambda(0:lm)
	real*8 mu(0:lm,0:lm),hz2,f(im+2,0:lm),fh(im+2,0:lm)
	real*8 left(im+2),righ(im+2),cent(im+2),hx,hx2
	real*8 fr1(0:lm),fr2(0:lm),fz1(im+2),fz2(im+2)
	real*8 cent1(im+2),left1(im+2),righ1(im+2)
	complex*16 wsave(3*lm),a(2*lm)
	integer i,k
	common/www/wsave

      print*,'in neumannsolver ' 
      stop
c        call zfft1d(a,2*lm,0,wsave)


      pi=3.14159265358979d0 

	hz2 = hz*hz
	hx2 = hx**2

	do i = 0,lm
         lambda(i) = 4d0/hz2*dsin(i*pi/2d0/lm)**2 !i!dcos(pi/lm*i)!1d0
	enddo

      do i = 1,im+2
         b(i,0)  = b(i,0)  - 2d0*fz1(i)*hx2
         b(i,lm) = b(i,lm) + 2d0*fz2(i)*hx2  
      enddo 

	do k = 0,lm
	   b(2,k)    = b(2,k)    - left(2)*fr1(k)*hx2
	   b(im+1,k) = b(im+1,k) - righ(im+1)*fr2(k)*hx2
	enddo 

	do i = 2,im+1
	   call directFT(b,fh,i)
	   left1(i) = left(i)*hx2
	   righ1(i) = righ(i)*hx2
      enddo

	do k = 0,lm
	   do i = 1,im+2
   	      cent1(i) = hx2*(cent(i) + lambda(k))
         enddo

         call crprogon(left1,cent1,righ1,u(1,k),fh(1,k))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      enddo 

	do i = 2,im+1
	   call inverseFT(u,f,i) 

    	   do k = 0,lm
	      u(i,k) = f(i,k)*(2d0/lm)  !/lm/im
	   enddo
      enddo
    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    	do k = 0,lm
	   u(1,k)    = fr1(k)
	   u(im+2,k) = fr2(k)
	enddo
    	
	end


      subroutine crprogon(a,b,c,x,f) 
      implicit none
      include 'four.par'
    
      real*8 x(im+2),f(im+2),a(im+2),c(im+2),b(im+2)
      real*8 al(im+2),be(im+2)
      integer i

      be(1)=f(1)/b(2) 
	al(1)=-c(1)/b(2)
	do i=2,im+2 
         al(i)  = -c(i)/(b(i-1)+al(i-1)*a(i)) 
         be(i)  = (f(i)-a(i)*be(i-1))/(b(i-1)+al(i-1)*a(i))
      enddo 

      x(im+2) = be(im+2)
      do i=im+1,2,-1 
         x(i)=al(i)*x(i+1)+be(i) 
      enddo 

      end

      subroutine crprogon1(a,b,c,x,f) 
      implicit none
      include 'four.par'
    
      real*8 x(im+2),f(im+2),a(im+2),c(im+2),b(im+2)
      real*8 al(im+2),be(im+2),hx2,t
      integer i

      t = x(im+2)

      be(2)=f(1)/b(2)
	al(2)=c(1)/b(2)
	do i=3,im+2 
         al(i)  = c(i-1)/(b(i-1)-a(i-1)*al(i-1)) 
         be(i)  = (f(i-1)+a(i-1)*be(i-1))/(b(i-1)-a(i-1)*al(i-1)) 
      enddo 

      x(im+2) = 0 ! be(im+2)
      do i=im+1,2,-1 
         x(i)=al(i+1)*x(i+1)+be(i+1) 
      enddo 

	x(im+2) = t

      end 

      subroutine directFT(b,fh,i)
      implicit none
      include 'four.par'
	real*8 h,b(im+2,0:lm),u(0:lm),u1(0:lm),pi,lambda(0:lm)
	real*8 mu(0:lm,0:lm),fh(im+2,0:lm),h2
	real*8 f(0:lm),rho(0:lm),fhw(0:lm)  
	integer i,k,k1,j

      pi=3.14159265358979d0  

      do k = 0,lm
	   fh(i,k) = 0

	      if((k.eq.0).or.(k.eq.lm)) then
		     rho(k) = 0.5d0*b(i,k)
            else
		     rho(k) = 1d0*b(i,k)
            endif
	enddo

      call wrapDirectFFTc(rho,fhw)

      do k = 0,lm
	   do j = 0,lm
c	      fh(i,k) = fh(i,k) + rho(j)*dcos(k*pi*j/lm)
	   enddo
	   fh(i,k) = fhw(k)
	enddo

      end

      subroutine inverseFT(b,fh,i)
      implicit none
      include 'four.par'
	real*8 h,b(im+2,0:lm),u(0:lm),u1(0:lm),pi,lambda(0:lm)
	real*8 mu(0:lm,0:lm),fh(im+2,0:lm),h2,fh1(im+2,0:lm)
	real*8 f(0:lm),rho(0:lm),fhw(0:lm) 
	integer i,k,k1,j

      pi=3.14159265358979d0  

      do k = 0,lm
	   fh(i,k) = 0

	      if((k.eq.0).or.(k.eq.lm)) then
		     rho(k) = 0.5d0*b(i,k)
            else
		     rho(k) = 1d0*b(i,k)
            endif
	enddo

      call wrapInverseFFTc(rho,fhw)

      do k = 0,lm
	   do j = 0,lm
c	      fh(i,k) = fh(i,k) + rho(j)*dcos(k*pi*j/lm)
	   enddo
	   fh(i,k) = fhw(k)
	enddo

      end


      subroutine wrapDirectFFTc(b,fh) 
      implicit none
      include 'four.par'

	real*8 b(0:lm),fh(0:lm)
	complex*16 a(2*lm)
	integer p,p1

      do p = 1,2*lm
	   a(p) = 0
	enddo

      do p=0,lm 
         a(p+1)=b(p)
      enddo 
      call fft(a,-1) 
      do p=1,lm+1 
         fh(p-1)=dreal(a(p)) 
      enddo 
	end

      subroutine wrapInverseFFTc(b,fh) 
      implicit none
      include 'four.par'

	real*8 b(0:lm),fh(0:lm)
	complex*16 a(2*lm)
	integer p,p1

      do p = 1,2*lm
	   a(p) = 0
	enddo

      do p=0,lm 
         a(p+1)=b(p)
      enddo 
      call fft(a,1) 
      do p=1,lm+1 
         fh(p-1)=2*lm*dreal(a(p)) 
      enddo 
	end

      subroutine fft(a,nsign)
	implicit none
	include 'four.par'
	complex*16 a(2*lm)
	integer nsign,k,j
	real*8 fh(0:lm),pi
	complex*16 wsave(3*lm)
	common/www/wsave

        pi=3.14159265358979d0   

	call fftc(a,nl,nsign,2*lm)

c        call zfft1d(a,2*lm,nsign,wsave)



	end
