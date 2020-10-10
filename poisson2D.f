c     NOTA BENE 'k' is absolute wavenumber, lk is PE-based
      subroutine poisson2D(k,lk) 
      implicit none 
c*************************************************
      integer puasfile,n_p,procs,parakm,im3,dlev
      integer rank,size 
c*************************************************
      	      
      integer im,im2,lm,n,i,k,l,km1,l2,lit 
      integer lit0 
      include 'pg2.par' 
      include 'para_pg2.par' 

      real*8 rr(2*im1)
      real*8 ph(2*im1,2*parakm,lm1),rh(im1,km),si2(km) 
      real*8 har(2*im1,lm1),si2seq(km),x(2*im1,lm1)
      real*8 pr(2*im1),c4dens(2*im1)
      real*8 cent(2*im1),righ(2*im1),left(2*im1),c3omit 
      real*8 s1,s2,s5,s6,eps,hr,hf,hz,tau,omit,pi,eps2,eps3 
      real*8 c1,c3,c4,ams0,hr2,hf2,eps1
      integer*4 ie(km),lep(2*parakm) 
      integer lk,prev,next,crank,harprocs,steady,group(procs),k0,i0,nt
      integer l0
      real*8 tol
      real*8 har1(2*im1,lm1),har2(2*im1,lm1),tol1,tol2
      !DEBUG ONLY: Eliminate!!
      real*8 pt1,ctex,tmp(2*im1) 
      
c     SWEEPING COEFS 
      real*8 c0pr(2*im1-1),apr(2*im1-1),bpr(2*im1-1),
     = al(2*im1-1),spr(2*im1-1),t
      
      common/compr/al,spr,apr,bpr,c0pr
      common/a/tau,hr,hf,hz,rr 
      common/g1/ph 
      common/pf/pr,c1,c3 
      common/si/si2 
      common/ts/ie,lep 
      common/x2/hr2,hf2,ams0,eps,omit 
      common/rhdens/rh
      common/harmonic/har
      common/grp/prev,next,crank,harprocs,group
      common/seq/si2seq
      common/aitken/har1,har2,tol1,tol2,steady
      common/ntc/nt
c**********************************************************       
      common/para1/rank,size

      print*, 'harmonic: global ',k,' local ',lk
      !stop

!      pt1 = amicro()
	  
      ctex = 0d0
c**********************************************************       
 
      eps1=eps/4d0 
      pi=3.14159265358979d0 
      lit=0 
      im2=im1+1 
      im=im1-1 
      km1=km2-1 
      lm=lm1-1 
      c1=1d0/hr**2 
      c3=1d0/hz**2 
      c4=4d0*pi/hz
      
      omit=1.25d0+0.75d0*ie(k)/dsqrt(8100d0+(ie(k))**2d0) 
      eps2=eps/dsqrt(si2(k)+1d0) 
      eps3=eps2/4d0


      !open(91,file='harm.txt',form='formatted')
c     ! 74  format(2i5,e25.15)     
c     PREPARING CURRENT HARMONIC
      do l=1,lm1 
          do i=2,2*im1 
	      har(i,l) = ph(i,lk,l)
              !write(91,74) i,l,har(i,l)
	  enddo
          har(1,l) = 0.0d0
      enddo
      !stop

      if(k.eq.-1) goto 2557
 
      
c     PREPARING COEFS FOR ITERATION
      do i=2,2*im1 
         cent(i)=(1.0d0-omit)*(2.0d0*(c1+c3)+ si2(k)/rr(i)**2) 
	 righ(i)=(1.0d0-omit)*c1*(i-1.0d0)/(i-1.5d0)
	 left(i)=(1.0d0-omit)*c1*(i-2.0d0)/(i-1.5d0)
      enddo
      c3omit = c3*omit
	 
c     COMPUTING SWEEPING COEFS       
      c0pr(2)=2.0d0*(c1+c3)+si2(k)/rr(2)**2
      al(2)  =2.0d0*c1/c0pr(2) 
      do i=3,2*im1-1                   
         c0pr(i)= 2.0d0*(c1+c3)+si2(k)/rr(i)**2
         apr(i) = (i-2.0d0)*c1/(i-1.5d0) 
         bpr(i) = (i-1.0d0)*c1/(i-1.5d0) 
         spr(i) = c0pr(i)-al(i-1)*apr(i) 
         al(i)  = bpr(i)/spr(i) 
      enddo 

c     PREPARING DENSITY ARRAY	  
      do i=2,2*im1 
        if(i.le.im1) then 
	   ! k IS AN ABSOLUTE WAVENUMBER
           c4dens(i)=c4*omit*rh(i,k) 
        else 
           c4dens(i)=0.0d0 
        endif 
      enddo
c##################################### for BCG testing
c 596  format(2i10)
c 597  format(i10,e25.16)
c598  format(i10,3e25.16)
c      
c      open(78,file='poisson.txt',form='formatted')
c      
c      write(78,596) 2*im,lm
c      
c      write(78,597) -1,c3
c      do i = 2,2*im1-1
c         write(78,597) i - 2,-c3*har(i,lm1) 
c      enddo  
c          
c      do i = 1,lm
c         write(78,597) i, -righ(2*im1 - 1)*har(2*im1,i)/(1.0d0-omit) 
c      enddo      
c
c      do i = 2,2*im1-1
c         write(78,597) i - 2,c4dens(i)/omit 
c      enddo      

c      do i = 2,2*im1-1
c         write(78,598) i - 2,left(i)/(1.0d0-omit),-cent(i)/(1.0d0-omit),
c     =	                 righ(i)/(1.0d0-omit) 
c      enddo      
c      close(78)
c      stop
c##################################### 

      lit0=0 
      n=0 
      s5=1d0 

   10 n=n+1 
      s1=0d0 
      s6=0d0 
c     GROUND LAYER       
          do 121 i=2,2*im+1 
	        pr(i)=2.0d0*c3omit*har(i,2)-c4dens(i)- 
     =  	      (righ(i)*har(i+1,1)+left(i)*har(i-1,1) 
     =                -cent(i)*har(i,1)) 
  121     continue 
  
    	  call progon(1,s1,s6) 
          l2=1 
          if(s6.lt.eps3) goto 17 
      print*,'after ground layer'
      !stop

c     HIGHER LAYERS 
      do 15 l=2,lm 
    		s6=0. 
    		do 141 i=2,2*im+1 
        	    pr(i)=c3omit*(har(i,l+1)+har(i,l-1))- 
     =                    (righ(i)*har(i+1,l)+left(i)*har(i-1,l) 
     =                     -cent(i)*har(i,l)) 
  141 	         continue 
    	         call progon(l,s1,s6) 
    	         l2=l 
                 
      if(s6.lt.eps3) goto 17 
	
   15 continue 
   17 continue

      print*,'iteration loop'
      !stop

      s2=abs(s5-s1)+1d-8 
      s2=s1*s5/s2 
      s5=s1 
      lit0=lit0+l2 
  101 format(i6,e25.15) 
      
      tol = s1
      
      if(n.lt.10) goto 10 

 1788 format('** ',i10,2e25.15)      
      if(tol.gt.eps3) goto 10

      if(dlev.ge.6) write(65,101) n,tol

c########################## CONTROL TRACING - STRONGLY NEEDED (09.06.03)
      n = n + 1
      s1 = 0.
      s6 = 0.
      
c     GROUND LAYER       
          do 221 i=2,2*im+1 
	        pr(i)=2d0*c3omit*har(i,2)-c4dens(i)- 
     =  	      (righ(i)*har(i+1,1)+left(i)*har(i-1,1) 
     =                -cent(i)*har(i,1)) 
  221     continue 
  
    	  call progon(1,s1,s6) 

c     HIGHER LAYERS 
      do 215 l=2,lm 
    		s6=0d0 
    		do 241 i=2,2*im+1 
        	    pr(i)=c3omit*(har(i,l+1)+har(i,l-1))- 
     =                (righ(i)*har(i+1,l)+left(i)*har(i-1,l) 
     =                -cent(i)*har(i,l)) 
  241 	    	continue 
    		call progon(l,s1,s6) 
	
  215 continue 

        tol = s1
 
      lit0=lit0+lm 
      if(tol.gt.eps2) goto 10 

c#######################################
      
      lit=lit+lit0 
      n=isign(1,lep(lk)) 
      if(lit0.ge.n*lep(lk)) then 
         n=-n 
      endif 

      ie(k)=ie(k)+n 
      lep(lk)=abs(n)*lit0 

    

      print*, 'harmonic ',k,' local ',lk,' done ',lep(lk)
      !stop
      k = -1

c##########################################
 2557 if(k.eq.-1) then
 
         do i=1,2*im1 
            cent(i)=(2.0d0*c1+ si2(k)/rr(i)**2) 
   	    righ(i)=-c1*(i-1.0d0)/(i-1.5d0)
	    left(i)=-c1*(i-2.0d0)/(i-1.5d0)
         enddo

         do i=2,2*im1 
            if(i.le.im1) then 
	   ! k IS AN ABSOLUTE WAVENUMBER
               c4dens(i)=2d0*pi/hz*rh(i,k) 
            else 
               c4dens(i)=0.0d0 
            endif 
         enddo
      ! endif

c      if(k.eq.1) then
         do l = 1,lm1
	      x(1,l) = har(1,l)
	      x(2*im1,l) = har(2*im1,l)
	   enddo  
         do i = 1,2*im1
	      x(i,lm1) = har(i,lm1)
	   enddo  
         print*,'neuwrap'
           !stop
	   call neuwrap(hr,hz,left,cent,righ,c4dens,x)

         do l = 1,lm1
	      x(1,l) = x(2,l)
	   enddo  

        t = 0
	do i = 1,im+1
c	   l = 1
    	   do l = 1,lm1
            if(t.lt.dabs(x(i,l)-har(i,l))) then
	         t = dabs(x(i,l)-har(i,l))
               i0 = i
	         l0 =l
	      endif
	      tmp(i) = x(i,1) - har(i,1)
	      har(i,l) = x(i,l)
	   enddo
        enddo
	print*,'nt,k,lk ',nt,k,lk
      endif
c##########################################

c     CURRENT HARMONIC DONE
      do i=1,2*im1 
          do l=1,lm1 
	      ph(i,lk,l) = har(i,l) 
	  enddo
c	  if((i.le.15).and.(k.eq.1)) write(65,350) i,har(i,1)
      enddo

      
      return 
      end 
      


c--------------------------------------------------------- 
      subroutine progon(l,s1,s6) 
      implicit none
c*************************************************
      integer puasfile,n_p,procs,parakm,im3
      integer rank,size,dlev 
c*************************************************
 
      include 'pg2.par' 
      include 'para_pg2.par' 
      
      real*8 har(2*im1,lm1)
      real*8 pr(2*im1),al(2*im1-1),be(2*im1-1)
      real*8 c1,c3,s1,s6,s2,s3
      integer i,l

      real*8 c0pr(2*im1-1),apr(2*im1-1),bpr(2*im1-1),
     = spr(2*im1-1)
      common/compr/al,spr,apr,bpr,c0pr
      common/para1/rank,size

      common/pf/pr,c1,c3 
      common/harmonic/har

      be(2)=pr(2)/c0pr(2) 
      do i=3,2*im1-1 
         be(i)=(pr(i)+apr(i)*be(i-1))/spr(i) 
      enddo 

      do i=2*im1-1,2,-1 
         s2=al(i)*har(i+1,l)+be(i) 
         s3=dabs(s2-har(i,l)) 
         s1=max(s1,s3) 
         s6=max(s6,s3) 
         har(i,l)=s2 
      enddo 
      har(1,l)=har(2,l) 
      return 
      end

