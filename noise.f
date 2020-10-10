      subroutine noise1 
      implicit none
      integer puasfile,n_p,procs,parakm,im3
      integer dlev
      include 'pg2.par' 
      include 'para_pg2.par' 

      integer im2,im,km1,lm,i,k,j,jm,nt,ml,n7,n8
      real*8 x,y,s1,s2,hr,hf,hz,qr1,tau,qr2
      real*8 f(100) 
      real*8 ro(im1+1,km2),rr(2*im1),ns(im1+1,km2) 
      real*8 xr(n_p),xf(n_p),vr(n_p),vf(n_p),qr(n_p) 
      common/a/tau,hr,hf,hz,rr 
      common/a1/f,ml,n7,n8 
      common/ntc/nt
      common/d1/ro 
      common/e1/jm,xr,xf,vr,vf
      common/e2/qr1,qr
      common/n/ns

      im2=im1+1 
      im=im1-1 
      km1=km2-1 
      lm=lm1-1 
      do 1 i=1,im2 
      do 1 k=1,km2 
      ns(i,k)=0d0 
    1 continue   
      do 2 j=1,jm 
      x=xr(j) 
      y=xf(j)
      qr2=qr(j)
      s1=x/hr+1.5d0 
      i=s1 
      s1=s1-i 
      s2=y/hf+1.5d0 
      k=s2 
      s2=s2-k 
          ns(i,k)=max(ns(i,k),qr2*(1d0-s1)*(1d0-s2))
          ns(i,k+1)=max(ns(i,k+1),qr2*(1d0-s1)*s2)
          ns(i+1,k)=max(ns(i+1,k),qr2*s1*(1d0-s2))
          ns(i+1,k+1)=max(ns(i+1,k+1),qr2*s1*s2)
    2 continue 
      do 3 i=1,im2 
      ns(i,2)=max(ns(i,2),ns(i,km2))
      ns(i,km1)=max(ns(i,km1),ns(i,1))
      ns(i,km2)=ns(i,2) 
      ns(i,1)=ns(i,km1) 
    3 continue 
      s1=(6d0*im+1d0)/(6d0*im-1d0) 
!      do 4 k=1,km2 
!      ro(2,k)=ro(2,k)-ro(1,k) 
!      ro(im1,k)=ro(im1,k)+s1*ro(im2,k) 
!    4 continue   
      do 5 k=1,km2 
      do 5 i=2,im1 
      if (ro(i,k).gt.0d0) ns(i,k)=ns(i,k)/ro(i,k) 
    5 continue 
!      do 6 k=1,km2 
!      ro(1,k)=ro(2,k) 
!      ro(im2,k)=ro(im1,k) 
!    6 continue 
      return 
      end 

      subroutine rad_char(pp,av_p,min_p,max_p)
      implicit none
      integer puasfile,n_p,procs,parakm,im3
      integer dlev
      include 'pg2.par' 
      include 'para_pg2.par' 

      integer im2,im,km1,lm,i,k,j,jm,nt,ml,n7,n8
      real*8 pp(im1+1,km2),av_p(im1),min_p(im1),max_p(im1)

      do i=1,im1
        av_p(i)=pp(i,1)
        min_p(i)=pp(i,1)
        max_p(i)=pp(i,1)
      do k=2, km2
       av_p(i)=av_p(i)+pp(i,k)
       min_p(i)=min(min_p(i),pp(i,k))
       max_p(i)=max(max_p(i),pp(i,k))
      enddo
       av_p(i)=av_p(i)/km2
      enddo
      end