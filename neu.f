c     hx,hz          - radial and "z" grid steps
c     left,cent,righ - radial stencil
c     fr1,fr2        - radial (Dirichlet) boundary conditions
c     fz1            - lower "z" (Neumann) boundary condition
c     fzul,fzuu      - two arrays for upper "z" (Neumann) boundary condition
c     b              - right-hand side
c     u              - solution 

c      subroutine neuwrap(hx,hz,left,cent,righ,
c     =                       	 fr1,fr2,fz1,fzul,fzuu,b,u1) 
      
      subroutine neuwrap(hx,hz,left,cent,righ,fz1,u1) 
 	implicit none
      include 'four.par'
	real*8 hz,b(im+2,0:lm),u1(im+2,0:lm),pi,lambda(0:lm),u(im+2,0:lm)
	real*8 mu(0:lm,0:lm),hz2,f(im+2,0:lm),fh(im+2,0:lm)
	real*8 left(im+2),righ(im+2),cent(im+2),hx,hx2
	real*8 fr1(0:lm),fr2(0:lm),fz1(im+2),fz2(im+2)
	real*8 fzul(im+2),fzuu(im+2)
	real*8 cent1(im+2),left1(im+2),righ1(im+2),t,fig2m1(im+2)
	real*8 bfr1(0:lm),bfr2(0:lm),bfz1(im+2),bfz2(im+2)
	integer i,k,i0,k0,nt

      common/g2m1/fig2m1
      common/ntc/nt
	common/bounds/bfz1,bfz2,bfr1,bfr2

      print*,'in neuwrap'
      !stop
      do i = 1,im+2
	   do k = 0,lm
	      b(i,k) = 0d0
	   enddo 
	enddo

	hx2 = hx*hx
	hz2 = hz*hz

      if(nt.eq.0) then
         do i = 1,im+2
c           fz1(i)  = (u1(i,1)-u1(i,0))/hz2 
            fz2(i)  = (u1(i,lm)-fig2m1(i))/hz2
	      bfz2(i) = fz2(i)
         enddo 

   	   do k = 0,lm
	      fr1(k) = u1(1,k)       !/hx2
	      fr2(k) = u1(im+2,k) !/hx2
	      bfr1(k) = fr1(k)       !/hx2
	      bfr2(k) = fr2(k) !/hx2
	      do i = 1,im+2
	         u1(i,k) = 0
	      enddo
	   enddo 
      else
         do i = 1,im+2
	      fz2(i)  = bfz2(i)
         enddo 

   	   do k = 0,lm
	      fr1(k) = bfr1(k)       !/hx2
	      fr2(k) = bfr2(k) !/hx2
	   enddo 
 
      endif

         do k = 0,lm
	      do i = 1,im+2
	         u1(i,k) = 0
	      enddo
	   enddo 

      print*,'b neumannsolver '
      !stop
      call neumannsolver(hx,hz,left,cent,righ,
     =                       	 fr1,fr2,fz1,fz2,b,u1) 

c      t = 0
c	do i = 2,im+1
c    	   do k = 0,lm-1
c         k = 0
c	      if(t.lt.dabs(u(i,k)-u1(i,k))) then
c	         t = dabs(u(i,k)-u1(i,k))
c               i0 = i
c	         k0 =k
c	      endif
c	   enddo
c      enddo

c	print*,t,t/dabs(u(i0,k0)-u(i0,k0))
c	stop

	end

