!  perevodit dvumerniy nomer yacheyki v odnomerney
!  postrochno
      subroutine dnumb(i,k,l)
      
      integer i,k,l
      include 'pg2.par'
      
      l=(k-1)*im1+i
      end

!  numb(i) soderzhit nomer posledney chastitsi v yacheyke
!  protsedura obnovlyaet massiv posle peredvizheniya i sortirovki     
      subroutine numbering(l1)
      implicit none
      integer puasfile,n_p,procs,parakm,im3
      integer rank,size,dlev
      integer i,j,l,l1,jm
      include 'pg2.par'
      include 'para_pg2.par'
      
      real*8 xr(n_p),xf(n_p),vr(n_p),vf(n_p),qr(n_p),qr1
      integer numcel(n_p),numb(0:im1*km2)
      common/num/numcel,numb            
      common/e1/jm,xr,xf,vr,vf
      common/e2/qr1,qr
      
!      print *,'numbering hello',jm
      l=im1*km2
      if (numb(0).ne.0) numb(0)=0
      do i=1,l
       j=numb(i-1)+1
!       write(62, *)'++',i,j,numcel(j)
!       if (j.lt.numb(i-1)) j=numb(i-1)+1
       if (j.gt.jm) then
         numb(i)=jm
       else
!        if (numcel(j).eq.i)then
        
         do while ((numcel(j).eq.i).and.(j.lt.jm))
          j=j+1
          enddo
!         write(62, *)'++',i,j,numcel(j) 
         if (j.eq.jm) then
           numb(i)=jm
         else
           numb(i)=j-1
         endif
       endif
!      print *,i,j
!          write(62, *),i,numb(i),numcel(numb(i))
      enddo
!          write(62, *),'=========================='
      end
      
      subroutine patincel(pic,x,y)
      implicit none
       integer puasfile,n_p,procs,parakm,im3
      integer rank,size,dlev
      
      real*8 x,y,tau,hr,hf,hz,rr,dr
      integer r,f,pic
      integer i,j,l
      include 'pg2.par'
      include 'para_pg2.par'
      
      integer numcel(n_p),numb(0:im1*km2)
      common/num/numcel,numb
      common/a/tau,hr,hf,hz,rr
      
      dr=x/hr+1d0
      r=dr
      dr=y/hf+1d0
      f=dr
      call dnumb(r,f,pic)
!      write(61, *),hr,hf,x,y,'<',r,f,'>',pic
      end
      
      subroutine amip(nb1)
      implicit none
      integer puasfile,n_p,procs,parakm,im3
      integer rank,size,dlev
            
      integer i,k,jm
      include 'pg2.par'
      include 'para_pg2.par'
      
      real*8 xr(n_p),xf(n_p),vr(n_p),vf(n_p),qr(n_p),qr1
      real *8 avm(im1*km2),amm
      integer numcel(n_p),numb(0:im1*km2)
      integer nu,nb,nb1
      common/num/numcel,numb
      common/e1/jm,xr,xf,vr,vf
      common/e2/qr1,qr
      common/massa/avm
      do i=1,im1*km2
        amm=0D0
        nu=numb(i-1)+1
        nb=numb(i)
        if (nu.le.nb)then
          do k=nu,nb
           amm=amm+qr(k)
          enddo
          amm=amm/(nb-nu+1)
        else
        amm=0
        endif
          avm(i)=amm
      enddo    
      end
      
      subroutine adaptivemass(j1)
      implicit none
      integer puasfile,n_p,procs,parakm,im3
      integer rank,size,dlev
            
      integer i,j,l,j1,jm,i5
      integer k,i1,k1,nrec,k2
      real*8 r16,aa,aa1,r4,a25
      include 'pg2.par'
      include 'para_pg2.par'
      
      real*8 xr(n_p),xf(n_p),vr(n_p),vf(n_p),qr(n_p),qr1
      real*8 ro(im1+1,km2),avm(im1*km2),rr(2*im1),tau,hr,hf,hz
      integer numcel(n_p),numb(0:im1*km2),hrt
      common/num/numcel,numb
      common/e1/jm,xr,xf,vr,vf
      common/e2/qr1,qr
      common/massa/avm
      common/d1/ro
      common/a/tau,hr,hf,hz,rr
      

      end