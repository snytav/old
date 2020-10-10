      subroutine finfile(n7) 
      implicit none
      integer puasfile,n_p,procs,parakm,im3,dlev
      integer rank,size 
      integer zh,fsave,fmom,fcen,fgas,fen,ffall,pot3d,fwin,
     =        pard,parfour,parv,gasd,gasv,rotc,rdisp,pdisp,dynflag,
     =        rvel,rflow
      real*8 cx,cy,dx,dy

      include 'pg2.par' 
      include 'para_pg2.par' 
      include 'auxillary.par'
      
      integer im2,km1,lm,i,k,n7,jm,j,k1,i1,v,nb,l
      real*8 pi,hr,c1,qr1,rm,r,pfi,u1,v1,x,y,u,s1,s2,vx,vy,sqdiff,
     =a,hf,tau,hz,hx,hy,xmin,xmax,ymin,ymax,di1,dk1,eps,s
      real*8 rr(2*im1),qr2 
      real*8 sigma(im1+1,km2),rf(im1,km),ro(im1+1,km2) 
      real*8 ur(im1+1,km2),uf(im1+1,km2),pp(im1+1,km2),urt 
      real*8 xr(n_p),xf(n_p),vr(n_p),vf(n_p),qr(n_p) 
      real*8 ux(im3,im3),uy(im3,im3)
      complex*16 rocx(km) 
c************************************************      
      real*8 fi(im1+1,km2),ams,aim,aen
      real*8 velx(im3,im3),rocart(im3,im3),vely(im3,im3),hr3
	real*8 tt(im1+1,km2),ns(im1+1,km2)
      real*8 velphi(im1),jmfull(im1),cphi(im1),crad(im1)
      real*8 velrad(im1)
      real*8 hat(7),hatf(7),hrt
      integer numcel(n_p),numb(0:im1*km2)
      real*8 av_p(im1),min_p(im1),max_p(im1)
      common/num/numcel,numb
      
      common/para1/rank,size
      common/para3/rocart,velx,vely,jmfull
      common/window/nb
c************************************************      
      character*12 tr,tr1 
      common/a/tau,hr,hf,hz,rr 
      common/b/ur,uf,pp 
      common/d/sigma 
      common/d1/ro
      common/e1/jm,xr,xf,vr,vf
      common/e2/qr1,qr
      common/e/fi,ams,aen,aim 
	common/temp/tt
	common/n/ns
	

c*************************************************
             
      if((pard.eq.0).and.(parv.eq.0).and.(parfour.eq.0).and.
     =   (gasd.eq.0).and.(gasv.eq.0).and.
     =	 (rotc.eq.0).and.(rdisp.eq.0).and.(pdisp.eq.0).and.
     =	 (pot3d.eq.0)) return

      if(dlev.ge.4) write(65,*) 'finfile started'	
 3022 format(i10,40x)    
      
c*************************************************
      pi=3.14159265358979d0 
      im2=im1+1 
      km1=km2-1 
      lm=lm1-1 
      
c**************************************************
      hr3=2.*hr*(im1 - 1.)/(im3 - 1.)

      c1=1d0/(hr3*hr3)
       
      do k=1,im3 
         do i=1,im3 
            rocart(i,k)   = 0. 
	    velx(i,k) = 0.
	    vely(i,k) = 0.
         enddo 
      enddo 
      
      do i=1,im1
        velphi(i) = 0.
	velrad(i) = 0.
	jmfull(i) = 0.
	cphi(i)   = 0.
	crad(i)   = 0.
      enddo
      
      if(rank.eq.0) then
          write(tr1,'(i3.3)') n7 
          if(pard.eq.1) then
	      tr='dpar'//tr1 
	      tr=tr(1:7)//'.dat' 
              open(80,file=tr,form='unformatted') 
          endif
	  
          if(parv.eq.1) then
	      tr='vpar'//tr1 
	      tr=tr(1:7)//'.dat' 
              open(81,file=tr,form='unformatted') 
          endif
	  
          if(gasd.eq.1) then
	      tr='dgas'//tr1 
	      tr=tr(1:7)//'.dat' 
              open(83,file=tr,form='unformatted') 
          endif
	  
          if(gasv.eq.1) then
	      tr='vgas'//tr1 
	      tr=tr(1:7)//'.dat' 
              open(84,file=tr,form='unformatted') 
          endif
	  
          if(pot3d.eq.1) then
	      tr='phi3d'//tr1 
	      tr=tr(1:8)//'.dat' 
              open(86,file=tr,form='unformatted') 
          endif

          if(parfour.eq.1) then
	      tr='fpar'//tr1 
	      tr=tr(1:7)//'.dat' 
              open(87,file=tr,form='unformatted') 
          endif
          
          if(rotc.eq.1) then
	      tr='vphi'//tr1 
	      tr=tr(1:7)//'.dat' 
              open(73,file=tr,form='formatted') 
          endif
	  
          if(pdisp.eq.1) then
	      tr='cphi'//tr1 
	      tr=tr(1:7)//'.dat' 
              open(74,file=tr,form='formatted') 
          endif

          if(rdisp.eq.1) then
	       tr='crad'//tr1
	       tr=tr(1:7)//'.dat'
	       open(75,file=tr,form='formatted') 
	  endif
	    tr='nois'//tr1 
	    tr=tr(1:7)//'.dat' 
            open(76,file=tr,form='formatted')    
	        
	    tr='n_rd'//tr1 
	    tr=tr(1:7)//'.dat' 
            open(77,file=tr,form='formatted')    
            
            tr='dp_r'//tr1 
	    tr=tr(1:7)//'.dat' 
            open(78,file=tr,form='formatted')    
	    
	    
      endif 	  
c**************************************************      	  
      if(dlev.ge.4) write(65,*) 'files opened'	

      rm=(im1-1)*hr 

      do j=1,jm 
         r=xr(j) 
         pfi=xf(j) 
         u1=vr(j) 
         v1=vf(j) 
         qr2=qr(j)
         x=r*dcos(pfi)+rm 
         y=r*dsin(pfi)+rm 
         u=u1*cos(pfi)-v1*sin(pfi) 
         v=u1*sin(pfi)+v1*cos(pfi) 
         s1=x/hr3+1. 
         i=s1              ! left 
         s1=s1-i 
         s2=y/hr3+1. 
         k=s2              ! left 
         s2=s2-k 
         rocart(i,k)=rocart(i,k)+(1.-s1)*(1.-s2)*qr2 
         rocart(i,k+1)=rocart(i,k+1)+(1.-s1)*s2*qr2 
         rocart(i+1,k)=rocart(i+1,k)+s1*(1.-s2)*qr2
         rocart(i+1,k+1)=rocart(i+1,k+1)+s1*s2*qr2 

         vx=cos(pfi)*vr(j)-sin(pfi)*vf(j) 
         vy=sin(pfi)*vr(j)+cos(pfi)*vf(j) 

         !///////
         velx(i,k)=velx(i,k)+(1.-s1)*(1.-s2)*vx 
         velx(i,k+1)=velx(i,k+1)+(1.-s1)*s2*vx 
         velx(i+1,k)=velx(i+1,k)+s1*(1.-s2)*vx 
         velx(i+1,k+1)=velx(i+1,k+1)+s1*s2*vx 

         !///////
         vely(i,k)=vely(i,k)+(1.-s1)*(1.-s2)*vy 
         vely(i,k+1)=vely(i,k+1)+(1.-s1)*s2*vy 
         vely(i+1,k)=vely(i+1,k)+s1*(1.-s2)*vy 
         vely(i+1,k+1)=vely(i+1,k+1)+s1*s2*vy 
	     
         s1=xr(j)/hr+1. 
         i=s1              ! left 
         s1=s1-i 

	 velphi(i)     = velphi(i) + (1 - s1)*vf(j)
	 velphi(i+1)   = velphi(i + 1) + s1*vf(j)

	 velrad(i)     = velrad(i) + (1 - s1)*vr(j)
	 velrad(i+1)   = velrad(i + 1) + s1*vr(j)

c        full number of particles is to be gathered
c        at each node
	 jmfull(i)     = jmfull(i) + (1. - s1)
	 jmfull(i + 1) = jmfull(i + 1) + s1
      enddo 
      if(dlev.ge.4) write(65,*) 'distributions computed'	
      
c*************************************************
      call addcart

      call addmeanvalue(velphi,velrad)
      
c     averaging velocities at each node
      do i = 1,im1
        if(jmfull(i).gt.0) then
    	    velphi(i) = velphi(i)/jmfull(i)
	    velrad(i) = velrad(i)/jmfull(i)
	endif
      enddo 
      if(dlev.ge.4) write(65,*) '... and prepared'	

      do j=1,jm 
         s1=xr(j)/hr+1.0d0 
         i=s1              ! left 
         s1=s1-i 

         sqdiff = (vf(j) - velphi(i))**2
	 cphi(i)     = cphi(i) + (1.0d0 - s1)*sqdiff
	 cphi(i+1)   = cphi(i + 1) +   s1*sqdiff

         sqdiff = (vr(j) - velrad(i))**2
	 crad(i)     = crad(i) + (1.0d0 - s1)*sqdiff
	 crad(i+1)   = crad(i + 1) +   s1*sqdiff
      enddo 

      call addmeanvalue(cphi,crad)
      
c     averaging dispersion values
      do i = 1,im1
        if(jmfull(i).gt.0) then
    	    cphi(i) = cphi(i)/jmfull(i)
	    crad(i) = crad(i)/jmfull(i)
	endif
      enddo 
      if(dlev.ge.4) write(65,*) 'curves prepared'	

      eps = 1.0d-08
      
c     NORMING VELOCITY VALUES
      do k=1,im3 
         do i=1,im3 
	    if(rocart(i,k).lt.eps) then
	       s = 1.0d0
	    else
	       s = rocart(i,k)
	    endif      
	    velx(i,k) = velx(i,k)/s
	    vely(i,k) = vely(i,k)/s
         enddo 
      enddo 


c writing binary distributions

      if(rank.eq.0) then	
	  hat(1)=nb         ! Variant number
	  hat(2)=im3
	  hat(3)=im3
	  hat(4)=-rm
	  hat(5)=-rm
	  hat(6)=hr3
	  hat(7)=hr3
	  
          if(pard.eq.1) then
	    write(80) (hat(i),i=1,7) 
	    write(80) ((c1*rocart(i,k),k=1,im3),i=1,im3) 
            close(80) 
 	  endif
          if(parv.eq.1) then
	    write(81) (hat(i),i=1,7) 
	    write(81) ((velx(i,k),k=1,im3),i=1,im3) 
	    write(81) ((vely(i,k),k=1,im3),i=1,im3) 
            close(81) 
	  endif

  100     format(2f10.4,e12.4) 
  101     format(2f10.4,2e12.4) 
  102     format(f10.4,e12.4) 

          if(rotc.eq.1) then
            do i = 1,im1
		write(73,102) (i - 0.5)*hr,velphi(i)
	    enddo
	    close(73)
	  endif     

          if(pdisp.eq.1) then
            do i = 1,im1
        	write(74,102) (i - 0.5)*hr,cphi(i)
	    enddo
	    close(74)
	  endif   

          if(rdisp.eq.1) then
              do i = 1,im1
		write(75,102) (i - 0.5)*hr,crad(i)
	    enddo
	    close(75)
	  endif 
	    call noise1
	    do i = 1,im1+1
	    do k=1,km2
	       call  dnumb(i,k,l)
		write(76,*)i,k, ns(i,k),numb(l)-numb(l-1)
	    enddo 
	    enddo
	    close(76)
	    call rad_char(ns,av_p,min_p,max_p)
	    do i=1,im1
		write(77,*) i*hr, av_p(i),min_p(i),max_p(i)
	    enddo
	    close(77)
	    call rad_char(ro,av_p,min_p,max_p)
	    do i=1,im1
		write(78,*) i*hr, av_p(i),min_p(i),max_p(i)
	    enddo
	    close(78)
      endif
      if(dlev.ge.4) write(65,*) 'particles written'	
	  
c************************************* here the gas begins      	  
      do k=1,im3 
         do i=1,im3 
            rocart(i,k)=0. 
            ux(i,k)=0. 
            uy(i,k)=0. 
         enddo 
      enddo 
      do k1=1,im3 
         y=(k1-0.5)*hr3-rm 
         do i1=1,im3 
            x=(i1-0.5)*hr3-rm 
            r=dsqrt(x*x+y*y) 
            if(r.ge.rm) then 
               a=0. 
               u=0. 
               v=0. 
            else 
               if((x.eq.0).and.(y.eq.0)) then
                 pfi = 0.0
               else 
               	 pfi=datan2(y,x) 
               endif 

               if(pfi.lt.0.) pfi=pfi+2.*pi 
               s1=r/hr+1.5 
               i=s1 
               s1=s1-i 
               s2=pfi/hf+1.5 
               k=s2 
               s2=s2-k 
               a=(1.-s1)*((1.-s2)*sigma(i,k)+s2*sigma(i,k+1))+ 
     =            s1*((1.-s2)*sigma(i+1,k)+s2*sigma(i+1,k+1)) 
               u1=(1.-s1)*((1.-s2)*ur(i,k)+s2*ur(i,k+1))+ 
     =            s1*((1.-s2)*ur(i+1,k)+s2*ur(i+1,k+1)) 
               v1=(1.-s1)*((1.-s2)*uf(i,k)+s2*uf(i,k+1))+ 
     =            s1*((1.-s2)*uf(i+1,k)+s2*uf(i+1,k+1)) 
               u=u1*cos(pfi)-v1*sin(pfi) 
               v=u1*sin(pfi)+v1*cos(pfi) 
            endif 
            rocart(i1,k1)=a 
            ux(i1,k1)=u 
            uy(i1,k1)=v 
         enddo 
      enddo
      if(dlev.ge.4) write(65,*) 'gas prepared'	

      if(rank.eq.0) then	

          if(gasd.eq.1) then
            write(83) (hat(i),i=1,7) 
      	    write(83) ((rocart(i,k),k=1,im3),i=1,im3) 
	    close(83)
	  endif

          if(gasv.eq.1) then
            write(84) (hat(i),i=1,7) 
	    write(84) ((velx(i,k),k=1,im3),i=1,im3) 
  	    write(84) ((vely(i,k),k=1,im3),i=1,im3) 
	    close(84)
	  endif
	  
      endif
      if(dlev.ge.4) write(65,*) 'gas written'	
      
      if((rvel.eq.1).and.(rank.eq.0)) then ! RADIAL VELOCITY vs RADIUS
         
	  tr='prvl'//tr1 
	  tr=tr(1:7)//'.dat' 
          open(80,file=tr,form='formatted') 
	  
  780     format(f15.8,f25.15)	  
	  do i = 2,im1
	     write(80,780) rr(i),velrad(i) 
	  enddo
	  close(80)
                  
	  tr='grvl'//tr1 
	  tr=tr(1:7)//'.dat' 
          open(80,file=tr,form='formatted') 
	  
	  do i = 2,im1
	     urt = 0.0d0
	     do k = 2,km1
	        urt = urt + ur(i,k) 
	     enddo
	     urt = urt / km
	     write(80,780) rr(i),urt 
	  enddo
	  close(80)
                  
      endif


      if((rank.eq.0).and.(dlev.eq.1)) then
  777    format(4e12.5)

	   tr='strm'//tr1 
	   tr=tr(1:7)//'.dat' 
         open(84,file=tr,form='formatted') 
	   do i = 2,im1
	      write(84,777) hr*(i-1.5d0),pp(i,2),sigma(i,2),tt(i,2)
	   enddo

	   close(84)
      endif 
	

      if((rflow.eq.1).and.(rank.eq.0)) then ! RADIAL FLOW vs RADIUS
         
	  tr='prfl'//tr1 
	  tr=tr(1:7)//'.dat' 
          open(80,file=tr,form='formatted') 
	  
	  do i = 2,im1
	     urt = 0.0d0
	     do k = 1,km1
	        urt = urt + ro(i,k)
	     enddo
	     urt = urt/km
	     write(80,780) rr(i),velrad(i)*urt 
	  enddo
	  close(80)
                  
	  tr='grfl'//tr1 
	  tr=tr(1:7)//'.dat' 
          open(80,file=tr,form='formatted') 
	  
	  do i = 2,im1
	     urt = 0.0d0
	     do k = 2,km1
	        urt = urt + ur(i,k)*sigma(i,k) 
	     enddo
	     urt = urt / km
	     write(80,780) rr(i),urt 
	  enddo
	  close(80)
                  
      endif

      if((parfour.eq.1).and.(rank.eq.0)) then
         hatf(1)=nb         ! Variant number
         hatf(2)=im1
         hatf(3)=km
         hatf(4)=0.0d0
	 hatf(5)=0.0d0
	 hatf(6)=hr
	 hatf(7)=hf
         
	 k1 = km/2 

         do i=2,im1 
            do k=1,km 
               rocx(k)=ro(i,k+1) 
            enddo 
            call fftc(rocx,nk,-1,km) 
            do k=1,k1+1 
               rf(i,k)=dreal(rocx(k)) 
            enddo 
            do k=2,k1 
               rf(i,k1+k)=dimag(rocx(k)) 
            enddo 
         enddo

         write(87) (hatf(i),i=1,7) 
      	 write(87) ((rf(i,k),k=1,km),i=1,im1) 
	 close(87)
         
      endif
c************************************* writing particles to the window
      if(fwin.ne.1) goto 5413

      hx = dx/im3
      hy = dy/im3
      if(dlev.ge.4) write(65,*) 'hx,hy,im3 ',hx,hy,im3	

      c1=1d0/(hx*hy)
       
      do k=1,im3 
         do i=1,im3 
            rocart(i,k)   = 0. 
	    velx(i,k) = 0.
	    vely(i,k) = 0.
         enddo 
      enddo 
      
      if(rank.eq.0) then
          write(tr1,'(i3.3)') n7 
	      tr='dpwin'//tr1 
	  tr=tr(1:8)//'.dat'
	   
	  if(pard.eq.1) then
              open(80,file=tr,form='unformatted') 
	  endif      

          write(tr1,'(i3.3)') n7 
	  tr='wparv'//tr1 
	  tr=tr(1:8)//'.dat' 
	  
	  if(parv.eq.1) then
              open(81,file=tr,form='unformatted') 
	  endif      
      endif 	  
c**************************************************      	  

      rm  = (im1 - 1.)*hr
      xmin = cx - dx/2. 
      ymin = cy - dy/2. 
      xmax = cx + dx/2. 
      ymax = cy + dy/2. 

      if(dlev.ge.4) then
         write(65,*) 'window files opened '
         write(65,*) 'hx,hy,xmin,ymin,xmax,ymax ',
     =	 hx,hy,xmin,ymin,xmax,ymax	
      endif	 



      hat(1)=nb         ! Variant number
      hat(2)=im3
      hat(3)=im3
      hat(4)=xmin
      hat(5)=ymin
      hat(6)=hx
      hat(7)=hy

      do j=1,jm 
         r=xr(j) 
         pfi=xf(j) 
         u1=vr(j) 
         v1=vf(j) 
         qr2=qr(j)
         x=r*cos(pfi) 
         y=r*sin(pfi)
	
	 if((x.gt.xmin).and.(x.lt.xmax).and.
     +	    (y.gt.ymin).and.(y.lt.ymax)) then  
         
            x = x - xmin
	    y = y - ymin

            u=u1*cos(pfi)-v1*sin(pfi) 
            v=u1*sin(pfi)+v1*cos(pfi) 
	    
	    s1=x/hx+1. 
            i=s1              ! left 
            s1=s1-i 
	    
	    s2=y/hy+1. 
            k=s2              ! left 
            s2=s2-k 
	    
	    if((i.lt.1).or.(i.ge.im3).or.
     +         (k.lt.1).or.(k.ge.im3)) then
		goto 497
     	    endif    
            rocart(i,k)=rocart(i,k)+(1.-s1)*(1.-s2)*qr2 
            rocart(i,k+1)=rocart(i,k+1)+(1.-s1)*s2*qr2 
            rocart(i+1,k)=rocart(i+1,k)+s1*(1.-s2)*qr2 
            rocart(i+1,k+1)=rocart(i+1,k+1)+s1*s2*qr2 

            vx=cos(pfi)*vr(j)-sin(pfi)*vf(j)    !# 
            vy=sin(pfi)*vr(j)+cos(pfi)*vf(j) 

            velx(i,k)=velx(i,k)+(1.-s1)*(1.-s2)*vx 
            velx(i,k+1)=velx(i,k+1)+(1.-s1)*s2*vx 
            velx(i+1,k)=velx(i+1,k)+s1*(1.-s2)*vx 
            velx(i+1,k+1)=velx(i+1,k+1)+s1*s2*vx 

           !///////
            vely(i,k)=vely(i,k)+(1.-s1)*(1.-s2)*vy 
            vely(i,k+1)=vely(i,k+1)+(1.-s1)*s2*vy 
            vely(i+1,k)=vely(i+1,k)+s1*(1.-s2)*vy 
            vely(i+1,k+1)=vely(i+1,k+1)+s1*s2*vy 

	 endif    
  497	 continue 
      enddo 

c*************************************************
      call addcart

      eps = 1.0d-08
c     NORMING VELOCITY VALUES
      do k=1,im3 
         do i=1,im3 
	    if(rocart(i,k).lt.eps) then
	       s = 1.0d0
	    else
	       s = rocart(i,k)
	    endif      
	    velx(i,k) = velx(i,k)/s
	    vely(i,k) = vely(i,k)/s
         enddo 
      enddo 

c writing binary distributions

      if(rank.eq.0) then	
	
          if(pard.eq.1) then  
	      write(80) (hat(i),i=1,7) 
	      write(80) ((c1*rocart(i,k),k=1,im3),i=1,im3) 
	      close(80)
	  endif
	  
	  if(parv.eq.1) then
	      write(81) (hat(i),i=1,7) 
	      write(81) ((velx(i,k),k=1,im3),i=1,im3) 
	      write(81) ((vely(i,k),k=1,im3),i=1,im3) 
              close(81)
	  endif      
      endif
      if(dlev.ge.4) write(65,*) 'wpart written'	

      if(rank.eq.0) then
          if(gasd.eq.1) then
              write(tr1,'(i3.3)') n7 
	      tr='dgwin'//tr1 
	      tr=tr(1:8)//'.dat' 
              open(83,file=tr,form='unformatted') 
	  endif      

          if(gasv.eq.1) then
              write(tr1,'(i3.3)') n7 
	      tr='wgasv'//tr1 
	      tr=tr(1:8)//'.dat' 
              open(84,file=tr,form='unformatted') 
	  endif      
      endif 	  

c**************************************************      	  
      do k=1,im3 
         do i=1,im3 
	    rocart(i,k)=0. 
            ux(i,k)=0. 
            uy(i,k)=0. 
         enddo 
      enddo

c*****writing gas fields to the window**************
       
      do k1=1,im3 
c    **** FIRST, GETTING GLOBAL COORDINATES OF A WINDOW POINT     
         dk1 = k1
         y=(dk1-0.5)*hy + ymin 
         do i1=1,im3 
            di1 = i1
            x=(di1-0.5)*hx + xmin 
            r=dsqrt(x**2+y**2) 
            if(r.ge.rm) then 
               a=0. 
               u=0. 
               v=0. 
            else 
               if((x.eq.0).and.(y.eq.0)) then
                 pfi = 0.0
               else 
               	 pfi=datan2(y,x) 
               endif 
               if(pfi.lt.0.) pfi=pfi+2.*pi 
               s1=r/hr+1.5 
               i=s1 
               s1=s1-i 
               s2=pfi/hf+1.5 
               k=s2 
               s2=s2-k 
               a=(1.-s1)*((1.-s2)*sigma(i,k)+s2*sigma(i,k+1))+ 
     =            s1*((1.-s2)*sigma(i+1,k)+s2*sigma(i+1,k+1)) 
               u1=(1.-s1)*((1.-s2)*ur(i,k)+s2*ur(i,k+1))+ 
     =            s1*((1.-s2)*ur(i+1,k)+s2*ur(i+1,k+1)) 
               v1=(1.-s1)*((1.-s2)*uf(i,k)+s2*uf(i,k+1))+ 
     =            s1*((1.-s2)*uf(i+1,k)+s2*uf(i+1,k+1)) 
               u=u1*cos(pfi)-v1*sin(pfi) 
               v=u1*sin(pfi)+v1*cos(pfi) 
            endif 
            rocart(i1,k1)=a 
            ux(i1,k1)=u 
            uy(i1,k1)=v 
         enddo 
      enddo
      if(dlev.ge.4) write(65,*) 'wgas prepared'	
      
      if(rank.eq.0) then	

          if(gasd.eq.1) then
            write(83) (hat(i),i=1,7) 
      	    write(83) ((rocart(i,k),k=1,im3),i=1,im3) 
	    close(83)
	  endif

          if(gasv.eq.1) then
            write(84) (hat(i),i=1,7) 
	    write(84) ((ux(i,k),k=1,im3),i=1,im3) 
  	    write(84) ((uy(i,k),k=1,im3),i=1,im3) 
	    close(84)
	  endif
      endif
      if(dlev.ge.4) write(65,*) 'wgas written'	

c****** writing potential to cartesian grid *****************
 5413 if(pot3d.eq.1) then
          call writephi(hr3,rm)
      endif
       
      return 
      end 

c---------------------------------------------------------      
      subroutine writephi(hr3,rm)
      implicit none 

      integer puasfile,procs,parakm,rank,size,im3
      integer n_p,im2,km1,dlev

      include 'pg2.par' 
      include 'para_pg2.par' 

      real*8 ph1(im1+1,km),ph2(im1+1,km),ph3(im1+1,km)
      real*8 hr3,rm,hz3,z1,phicart(im3,im3,5),t
      integer l,l1,l2,i,k,i1,k1	
      real*8 tau,hr,hf,hz,rr(2*im1),fi(im1+1,km2)
      real*8 x,y,r,pfi,s1,s2,pi,hat(10)
      complex*16 a(km)
      
      common/a/tau,hr,hf,hz,rr 
      common/para1/rank,size
      common/paraphi3d/ph1,ph3
      
      pi=3.14159265358979d0 
      
      im2 = im1 + 1
      km1 = km2 - 1

c     first, determining the output z step
      hz3 = ((lm1 - 1.)*hz)/(im3-1.)
      
      do l2=1,im3-1
	    
c       determining the height where we are now
        z1 = ((l2 - 1.)*hz3)/hz + 1.
	l1 = z1
	
c       obtaining interpolation factor 	
	z1 = z1 - l1 

c       gathering the lower and upper layers
c	call gatherLayers(l2) 

c       doing interpolation
        do i = 1,im1+1
	    do k = 1,km
		ph2(i,k) = (1. - z1)*ph1(i,k) + z1*ph3(i,k)
	    enddo
	enddo  
	
c       reverse Fourier transform	
	k1 = km/2
        do i=2,im2 
    	    do k=1,k1+1 
        	a(k)=dcmplx(ph2(i,k),0.d0) 
            enddo 
            do k=2,k1 
        	a(k)=a(k)+dcmplx(0.d0,ph2(i,k1+k)) 
            enddo 
         
	    do k=2,k1 
        	a(k1+k)=dconjg(a(k1+2-k)) 
            enddo 
         
	    call fftc(a,nk,1,km) 
            do k=1,km
		fi(i,k+1)=dreal(a(k)) 
            enddo 
         enddo 
      
c        fixing boundary values      
    	 do i=2,im2 
            fi(i,1)=fi(i,km1) 
            fi(i,km2)=fi(i,2) 
         enddo 
         do k=1,km2 
            fi(1,k)=fi(2,k) 
         enddo 
	
        do k=1,im3 
    	     do i=1,im3 
        	phicart(i,k,l2)=0. 
             enddo 
        enddo 
      
c        writing to cartesian grid      
        do k1=1,im3 
             y=(k1-0.5)*hr3-rm 
             do i1=1,im3 
                x=(i1-0.5)*hr3-rm 
                r=dsqrt(x*x+y*y) 
                if(r.ge.rm) then 
                   t=0. 
                else 
                   if((x.eq.0).and.(y.eq.0)) then
                      pfi = 0.0
                   else 
                      pfi=datan2(y,x) 
                   endif 
                   if(pfi.lt.0.) pfi=pfi+2.*pi 
                   s1=r/hr+1.5 
                   i=s1 
                   s1=s1-i 
                   s2=pfi/hf+1.5 
                   k=s2 
                   s2=s2-k 

                   t=(1.-s1)*((1.-s2)*fi(i,k)+s2*fi(i,k+1))+ 
     =                 s1*((1.-s2)*fi(i+1,k)+s2*fi(i+1,k+1)) 
                endif 
		phicart(i1,k1,l2) = t
             enddo 
        enddo 
      enddo     !vertical loop	

      if(rank.eq.0) then

        hat(1)=-10001         ! Variant number
        hat(2)=im3
	hat(3)=im3
	hat(4)=im3 
	hat(5)=-rm
	hat(6)=-rm
	hat(7)=0.
	hat(8)=hr3
	hat(9)=hr3
	hat(10)=hz3

	write(86) (hat(i),i=1,10) 
      	write(86) (((phicart(i,k,l),k=1,im3),i=1,im3),l=1,im3) 
	close(86)
     
      endif
c************************************************************
      
      end      
