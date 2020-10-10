c           of particle impulse(s3), total moment of impulse(s10)       
c    pgb -- gravity centre of gas (x,y)(s8,s9), gravity centre of  
c           particle (x,y)(s11,s12), total centre of gravity (x,y)(s13,s14) 
c    pgc -- kinetic gas energy(s6), inner energy of gas(s15),  
c           gas-field energy(s16) 
c    pgd -- particle kinetic energy(s4) and field energy(s17),  
c           total energy(s1) 
c     
c --------------------------------------------- 
      implicit none
      external amicro 
      integer puasfile,n_p,procs,parakm,im3
      integer dlev,zh
      integer rank,size 

      real*8  cx,cy,dx,dy,amicro
      integer fsave,fmom,fcen,fgas,fen,ffall,pot3d,fwin,
     =        pard,parfour,parv,gasd,gasv,rotc,rdisp,pdisp,dynflag,
     =        rvel,rflow

      include 'pg2.par' 
      include 'para_pg2.par' 
      include 'auxillary.par'
      real*8 fi(im1+1,km2),sigma(im1+1,km2), 
     =ro(im1+1,km2),rr(2*im1),f(100),tau,omit,hr,hf,tx,hz,zm,
     =rm,rd1,rd2,ams,aen,aim,hr2,hf2,ams0,eps
      real*8 ur(im1+1,km2),uf(im1+1,km2),si2(km),
     =pp(im1+1,km2),ur1(im1+1,km2),uf1(im1+1,km2) 
      character*4 st 
      character*12 tr 
      integer nt0,ml,nt1,nt2,m1,im2,km1,im,lm,m7,n8,nt,n7,m5 
      integer lit2,nb
      common/a/tau,hr,hf,hz,rr 
      common/a1/f,ml,n7,n8
      common/a2/zm,rd2,lit2,st 
      common/ntc/nt
      common/crm/rm,rd1
      common/b/ur,uf,pp 
      common/c/ur1,uf1 
      common/d/sigma 
      common/d1/ro 
      common/si/si2
      common/e/fi,ams,aen,aim 
      common/x2/hr2,hf2,ams0,eps,omit
c***************************************************
      common/para1/rank,size
      common/window/nb
c save all the data, flags for output: momentum, gravity centre movement,
c gas invariants, various energy, fall to the centre      
c      common/flags/fsave,fmom,fcen,fgas,fen,ffall,pot3d,fwin,
c     =      pard,parvx,parvy,gasd,gasvx,gasvy,rotc,rdisp,pdisp      

c     ZERO HARMONIC COMPUTATION FLAG
      ml = zh
              
      call initialize(ml)    
c***************************************************

      open(12,file='pgstart9.dat',form='formatted') 
      read(12,302) nt1 
      read(12,302) nt2 
      read(12,301) tau 
      read(12,301) omit
      close(12) 

      close(12)
 
  301 format(f10.4,40x)    
  302 format(i10,40x)    


      if(dlev.ge.1) write(65,*) 'read physical data'
      m1=12    
 9001 format(4x,'nt1,nt2=',2i4,4x,'ml,m1=',2i3,4x,'tau=',f11.4) 
 9003 format(4x,'im1=',i4,4x,'km2=',i4,4x,'lm1=',i4) 
      im2=im1+1 
      km1=km2-1 
      im=im1-1 
      lm=lm1-1 
      if(ml.gt.1) then
    	    open(89,file='check.txt',form='formatted') 
    	    read(89,1090) m7
 1090       format(i4.4)      
            close(89)
	    write(st,1090) m7
            nt0 = m7
      else 
         st='0000'
	 nt0 = 0
      endif 
      n8=0

      if(dlev.ge.2) write(65,*) 'preparing output files'      
      if(rank.eq.0) then
          tr='pg2'//st 
          tr=tr(1:7)//'.lst' 
          open(25,file=tr,form='formatted') 

          write(25,9001) nt1,nt2,ml,m1,tau 
          write(25,9003) im1,km2,lm1 
          write(25,*) 'hr,hf=',hr,hf 

          tr='pga'//st
	  tr=tr(1:7)//'.txt' 
          if(fmom.eq.1) open(31,file=tr,form='formatted') 
          
	  tr='pgb'//st 
          tr=tr(1:7)//'.txt' 
          if(fcen.eq.1) open(32,file=tr,form='formatted') 
          
	  tr='pgc'//st 
          tr=tr(1:7)//'.txt' 
          if(fgas.eq.1) open(33,file=tr,form='formatted') 
          
	  tr='pgd'//st 
          tr=tr(1:7)//'.txt' 
          if(fen.eq.1) open(34,file=tr,form='formatted') 

	  tr='ctr'//st 
          tr=tr(1:7)//'.txt' 
          if(ffall.eq.1) open(35,file=tr,form='formatted') 
          
	  write(25,9002) nt,tx 
	  
      endif
      if(dlev.ge.2) write(65,*) 'output files prepared'      

      call start(tx) 
   
      call printTime(nt,0)

      call energy(tx) 
             
c--------------------------------------------------------- 
      n8=n8+1 
      if(ml.gt.1) goto 1 
      call finfile(0) 
c      goto 4 
    1 do 3 m7=nt0+1,nt2+nt0
      n7=n7+1 
      do 2 m5=1,nt1 
      nt=nt+1 
      tx=tx+tau 
      if(dlev.ge.1) write(65,*) 
     =nt,' ====================================================='
 9002 format(/8x,'nt=',i8,1x,' tx=',f7.4) 

      call particle(m1) 

      call gaz(ams) 

c      call chem

      call reduce(nt,tx,ffall) 

      call energy(tx) 
      print*,'before puas7 ',nt
      call puas7(lit2) 

      call printTime(nt,0)

      if(rank.eq.0) then
        write(25,9002) nt,tx 
      endif

    2 continue   

      call finfile(n7)

      if(ml.ne.1.and.ml.ne.3) goto 3 

    3 continue 
    4 continue 
c---------------------------------------------------------       

      if(fsave.eq.1) call final(m7,tx) 

      close(25) 
      
      if(fmom.eq.1) close(31)
      if(fcen.eq.1) close(32)
      if(fgas.eq.1) close(33)
      if(fen.eq.1)  close(34)
      if(ffall.eq.1) close(35)

c***************************************************
      call finalize

      call printTime(nt,1)

c***************************************************
      
      end 

c***************************************************

c---------------------------------------------------------       
      subroutine start(tx) 
      implicit none
      integer puasfile,n_p,procs,parakm,im3
      integer rank,size,dlev
      real*8 cx,cy,dx,dy
      integer zh,fsave,fmom,fcen,fgas,fen,ffall,pot3d,fwin,
     =        pard,parfour,parv,gasd,gasv,rotc,rdisp,pdisp,dynflag,
     =        rvel,rflow
      
      include 'pg2.par' 
      include 'para_pg2.par' 
      include 'auxillary.par'

c***********************************************************      
      real*8 si2seq(km)
      real*8 surph(im1+1,km)
c***********************************************************      

      real*8 fi(im1+1,km2),fig1(lm1),fig2(2*im1),fig2m1(2*im1)  
      real*8 ph(2*im1,2*parakm,lm1),si2(km) 
      real*8 fr(im1+1,km2),ff(im1+1,km2) 
      real*8 xr(n_p),xf(n_p),vr(n_p),vf(n_p),qr(n_p) 
      real*8 ro(im1+1,km2),sigma(im1+1,km2),rr(2*im1),f(100) 
      real*8 ur(im1+1,km2),uf(im1+1,km2),pp(im1+1,km2), 
     =ur1(im1+1,km2),uf1(im1+1,km2),vvr(im1+1,km2),vvf(im1+1,km2) 
      real*8 rr1(im1+1),rr2(im1+1) 
      real*8 rm,rd1,tx,zm,rd2,ams0,hr,hf,hz,gamma,cg,tr,pi,
     =amp,amg,ams,aen,aim,eps,s1,tau,omit,vd,pp0,tr0,s,x,hr2,hf2,
     =r0,h1,s2,s3,r,a1,zz,s4,ro01,sigma0,omega,s8,sx,sy,s7,y,xf0,s9,
     =s10,er,us,vs,cs,sn,hz2,c1,s6,f9,t05dde,vx,vy,tau1,qr1
      integer*4 ie(km),lep(2*parakm)
      integer lit2,im,lm,im2,km1,ml,i,nb,nt,n7,jm,icen,k,l,j,k1,
     =nt1,nt2,icn,i2,nr,k2,j1,n8,i1
      real*8 g05cae,ellip,dabs,press,g05dde,tt(im1+1,km2),kappa,tt0
      complex*16 a(km) 
      character*4 st 
      character*15 tr2,tr1 

	real*8 cc(mm,im1+1,km2),cn(mm)
	
	integer mu(mm),m

      real*8 ttloc(im1+1,km2),
     &	ttloc_box

	real*8 qq(im1+1,km2),qq_box
     
	common/d11/cc

	common/d12/mu
	common/d13/ttloc,ttloc_box
	common/d14/qq,qq_box


      common/a/tau,hr,hf,hz,rr 
      common/a1/f,ml,n7,n8 
      common/ntc/nt
      common/a2/zm,rd2,lit2,st 
      common/crm/rm,rd1
      common/b/ur,uf,pp 
      common/c/ur1,uf1 
      common/d/sigma 
      common/d1/ro 
      common/d2/tr,vvr,vvf 
      common/e/fi,ams,aen,aim 
      common/e1/jm,xr,xf,vr,vf
      common/e2/qr1,qr
      common/g/fr,ff 
      common/g1/ph 
      common/i/icen 
      common/si/si2 
      common/ts/ie,lep 
      common/x2/hr2,hf2,ams0,eps,omit 
      common/x3/amp,amg 
 
c***********************************************************      
      integer harmonics(procs),up(procs),down(procs),group(procs)
      integer lodge(procs,2*parakm),workload(procs),layers(km)
      integer prev,next,crank,harprocs
      real*8  threshold
      common/lod/lodge
      common/wor/workload
      common/thd/threshold
      common/lay/layers
      common/har/harmonics
      common/harbounds/up,down
      common/grp/prev,next,crank,harprocs,group

      common/para1/rank,size
      common/para2/surph
      common/window/nb
      integer evn,m7

      common/g2m1/fig2m1
      common/x4/gamma,kappa
c***********************************************************      
      jm=0.9*n_p/procs
      if(dlev.ge.1) write(65,*) 'start procedure begins'

      lit2=0 
      im=im1-1 
      lm=lm1-1 
      im2=im1+1 
      km1=km2-1 
      if(ml.le.1) goto 4

      open(89,file='check.txt',form='formatted')
 1091 format(i4.4)
      read(89,1091) m7 
      close(89)

      evn = m7 - (m7/2)*2
      
      write(tr2,'(i2.2,i1.1)') rank,evn 
      tr1='ptot'//tr2(1:3)//'.dat' 

      open(10,file=tr1,form='unformatted') 
      read(10) (f(i),i=1,100) 
      tx=f(1) 
      nb=f(2)+0.1
      nt=f(3)+0.1 
      im=f(4)+0.1 
      lm=f(6)+0.1 
      rm=f(7) 
      zm=f(8) 
      rd1=f(9) 
      rd2=f(10) 
      ams0=f(11) 
      qr1=f(12) 
      hr=f(13) 
      hf=f(14) 
      hz=f(15) 
      n7=f(16)+0.1 
      gamma=f(17) 
      cg=f(18) 
      tr=f(19) 
      amp=f(20) 
      amg=f(21) 
      ams=f(22) 
      aen=f(23) 
      aim=f(24) 
      jm=f(25)+0.1 
      icen=f(26)+0.1 
      eps=f(27) 
c      cx = f(31)
c      cy = f(32)
c      dx = f(33)
c      dy = f(34)
      threshold = f(35)
      prev = f(36) + 0.1
      next = f(37) + 0.1
      crank = f(38) + 0.1
      harprocs = f(39) + 0.1

c     BALANCING DATA
      read(10) (harmonics(k),k=1,procs)
      read(10) ((lodge(i,k),i=1,procs),k=1,2*parakm)
      read(10) (workload(k),k=1,procs)
      read(10) (layers(k),k=1,km) 
      read(10) (up(k),k=1,procs)
      read(10) (down(k),k=1,procs)
      read(10) (group(k),k=1,procs)

      read(10) (ie(k),k=1,km) 
      read(10) (lep(k),k=1,2*parakm) 
      read(10) (rr(i),i=1,2*im1) 
      read(10) (si2(k),k=1,km) 
      read(10) (((ph(i,k,l),i=1,2*im1),k=1,2*parakm),l=1,lm1) 
      read(10) ((sigma(i,k),i=1,im2),k=1,km2) 
      read(10) ((pp(i,k),i=1,im2),k=1,km2) 
      read(10) ((ur(i,k),i=1,im2),k=1,km2) 
      read(10) ((uf(i,k),i=1,im2),k=1,km2) 
      read(10) (xr(j),j=1,jm) 
      read(10) (xf(j),j=1,jm) 
      read(10) (vr(j),j=1,jm) 
      read(10) (vf(j),j=1,jm) 
      close(10) 
c-------------------- восстановление потенциала       
      if(rank.eq.0) open(61,file='balance.txt',form='formatted')

      if(dlev.ge.2) write(65,*) 'read saved data'

c***********************************************      
      do i=1,im1+1
        do k=1,parakm
	    surph(i,k + rank*parakm) = ph(i,k,1)
	enddo
      enddo

      call gatherFi
c***********************************************      
      k1=km/2

      do i=2,im2 
         do k=1,k1+1 
            a(k)=dcmplx(surph(i,k),0.0d0) 
         enddo 
         do k=2,k1 
            a(k)=a(k)+dcmplx(0.0d0,surph(i,k1+k)) 
         enddo 
         do k=2,k1 
            a(k1+k)=dconjg(a(k1+2-k)) 
         enddo 
         call fftc(a,nk,1,km) 
         do k=1,km
	     fi(i,k+1)=dreal(a(k)) 
         enddo 
      enddo 
      
      do i=2,im2 
         fi(i,1)=fi(i,km1) 
         fi(i,km2)=fi(i,2) 
      enddo 
      do k=1,km2 
         fi(1,k)=fi(2,k) 
      enddo 

      do i=2,im1 
         s1=(ams0+ams)/((i-1.0d0)*hr)**2 
         do k=1,km2 
            fr(i,k)=(fi(i,k)-fi(i+1,k))/hr-s1 
	 enddo
      enddo	     

      do k=1,km2 
         fr(1,k)=fr(2,k) 
      enddo 
      do 20 i=1,im2   
      do 20 k=1,km1 
      ff(i,k)=(fi(i,k)-fi(i,k+1))/(hf*rr(i)) 
   20 continue   
      goto 14 
c==================================================       
    4 continue 

      open(11,file='conc04.sta',form='formatted')
      s=0.d0
      do m=1,mm
         read(11,301) cn(m)
         s=s+cn(m)
      enddo
      do m=1,mm
         cn(m)=cn(m)/s
      enddo
	
      open(13,file='mu01.sta',form='formatted')
	do m=1,mm
	 read(13,302) mu(m)
      enddo

      pi=3.14159265358979d0 

      open(12,file='pgstart9.dat',form='formatted') 
      read(12,302) nt1 
      read(12,302) nt2 
      read(12,301) tau 
      read(12,301) omit
      
      read(12,302) nb 
      read(12,301) rm 
      read(12,301) zm 
      read(12,301) rd1     ! радиус диска частиц 
      read(12,301) amp     ! масса частиц 
      read(12,301) vd      ! радиальная дисперсия частиц 
      read(12,301) rd2     ! радиус газового диска 
      read(12,301) amg     ! масса газа 
      read(12,301) tt0     ! давление в центре 
      read(12,301) gamma 
      read(12,301) kappa   ! коэффициент теплопроводности
      read(12,301) tr0     ! коэффициент трения 
      read(12,301) ams0     ! масса Солнца 
      read(12,302) nr 
      read(12,302) icn 
      read(12,301) eps 

c***************************************************
c      read(12,301) cx 
c      read(12,301) cy 
c      read(12,301) dx 
c      read(12,301) dy 
      
      close(12)
            
      call synchronize
     
c***************************************************

      icen=icn+1  
      do k=1,km 
         ie(k)=90.0d0*(4.0d0*omit-5.0d0)/
     =	       dsqrt(16.0d0*(omit-0.5d0)*(2.0d0-omit)) 
         if(k.lt.(2*parakm)) lep(k)=1 
      enddo 
      
      tr=tr0 
      ams=0.0d0 
      aen=0.0d0 
      aim=0.0d0 
      qr1=amp/(jm*procs) 
  301 format(f10.4,40x)    
  302 format(i10,40x) 
      s=0.0d0 
      x=1.0d0 
      do 27 i=1,nr 
      s=s+g05cae(x) 
   27 continue 

      hr=rm/im 
      hf=2.0d0*pi/km 
      hz=zm/lm 
      hr2=hr**2 
      hf2=hf**2 
      hz2=hz**2 
      tx=0.0d0 
      nt=0 
      n7=0 
c      n8=0 

      f(1)=tx 
      f(2)=nb 
      f(3)=nt 
      f(4)=im 
      f(5)=km 
      f(6)=lm 
      f(7)=rm 
      f(8)=zm 
      f(9)=rd1 
      f(10)=rd2 
      f(11)=ams0 
      f(12)=qr1 
      f(13)=hr 
      f(14)=hf 
      f(15)=hz 
      f(17)=gamma 
      f(19)=tr 
      f(20)=amp 
      f(21)=amg 
      f(22)=ams 
      f(23)=aen 
      f(24)=aim 
      f(25)=jm 
      f(26)=icen 
      f(27)=eps 
      do 5 i=1,2*im1 
         rr(i)=hr*(i-1.5d0) 
    5 continue 
c ============= гран.условие потенциала 1 ===================== 
      if(dlev.ge.2) write(65,*) 'computing boundary conditions'

      c1=-4.0d0/(pi*rm*rm) 
      r0=2.0d0*rm+hr/2.0d0 
      do 23 l=1,lm1 
      s=(hz*(l-1.0d0))**2 
      h1=rm/8.0d0 
      i2=8 
      s3=0.0d0 
   21 s2=0.0d0 
      r=-h1/2d0 
      do 22 i=1,i2 
      r=r+h1 
      s1=(r0+r)**2+s 
      a1=4.*r0*r/s1 
      s2=s2+r*ellip(a1)/dsqrt(s1)
   22 continue 
      s2=s2*h1 
      s1=dabs(s2-s3) 
      s3=s2 
      h1=h1/2.0d0 
      i2=2*i2 
      if(s1.gt.eps) goto 21 
      fig1(l)=c1*s2
   23 continue 
      if(dlev.ge.2) write(65,*) 'computed R boundary conditions'
   
c ============= гран.условие потенциала 2 ===================== 
      c1=-4.0d0/(pi*rm*rm) 
      zz=(zm-hz)**2 
      do 126 i=2,2*im1 
      r0=hr*(i-1.5d0) 
      h1=rm/8.0d0 
      k2=8 
      s3=0.0d0 
  124 s2=0.0d0 
      r=-h1/2.0d0 
      do 125 k=1,k2 
      r=r+h1 
      s1=(r0+r)**2+zz 
      a1=4.0d0*r0*r/s1 
      s2=s2+r*ellip(a1)/dsqrt(s1)
  125 continue 
      s2=s2*h1 
      s1=dabs(s2-s3) 
      s3=s2 
      h1=h1/2.0d0
      k2=2*k2 
      if(s1.gt.eps) goto 124 
      fig2(i)=c1*s2 
  126 continue 
      fig2(1)=fig2(2) 
      do i = 1,2*im1
	   fig2m1(i) = fig2(i) ! EXPERIMENTAL TO SET NEUMANN BOUNDARY CONDITIONS
	enddo  

      if(dlev.ge.2) write(65,*) 'computed Z boundary conditions'
c==============================================================
      c1=-4.0d0/(pi*rm*rm) 
      zz=zm*zm 
      do 26 i=2,2*im1 
      r0=hr*(i-1.5d0) 
      h1=rm/8.0d0 
      k2=8 
      s3=0.0d0 
   24 s2=0.0d0 
      r=-h1/2.0d0 
      do 25 k=1,k2 
      r=r+h1 
      s1=(r0+r)**2+zz 
      a1=4.0d0*r0*r/s1 
      s2=s2+r*ellip(a1)/dsqrt(s1) 
   25 continue 
      s2=s2*h1 
      s1=dabs(s2-s3) 
      s3=s2 
      h1=h1/2.0d0
      k2=2*k2 
      if(s1.gt.eps) goto 24 
      fig2(i)=c1*s2 
   26 continue 
      fig2(1)=fig2(2) 
      if(dlev.ge.2) write(65,*) 'computed Z boundary conditions'

c==============================================================      
      do i=1,im2 
         s1=1.0d0-(rr(i)/rd1)**2 
         if(s1.gt.0.) then 
            rr1(i)=dsqrt(s1) 
         else 
            rr1(i)=0.0d0 
         endif 
         s1=1d0-(rr(i)/rd2)**2 
         if(s1.gt.0.) then 
            rr2(i)=dsqrt(s1) 
         else 
            rr2(i)=0.0d0
         endif 
         do k=1,km2 
            ro(i,k)=rr1(i) 
            sigma(i,k)=rr2(i) 
	 enddo
      enddo 	      

      s1=0.0d0 
      
      if(dlev.ge.1) write(65,*) 'preparing distribution'
      
      s3=0. 
      do 2 i=2,im1 
         s2=0.0d0 
         s4=0.0d0 
         do 1 k=2,km1 
            s2=s2+ro(i,k) 
            s4=s4+sigma(i,k) 
    1    continue 
         s1=s1+s2*rr(i) 
         s3=s3+s4*rr(i) 
    2 continue 
    
      s1=hr*hf*s1 
      s3=hr*hf*s3 
c ----- коррекция плотности частиц и газа 
c ----- и вычисление давления 
      ro01=amp/s1 
      sigma0=amg/s3 
 
      if(rank.eq.0) then
         write(25,*) 'particle density in centre =',ro01 
         write(25,*) '     gas density in centre =',sigma0 
      endif

      do i=1,im2 
         do k=1,km2 
            ro(i,k)=ro01*ro(i,k) 
            sigma(i,k)=sigma0*sigma(i,k) 
            do m=1,mm
               cc(m,i,k)=cn(m)
   61       enddo
	 enddo
      enddo	     

      do 9 i=1,im2 
      do 9 k=1,km2 
      pp(i,k)=press(sigma(i,k)) 
    9 continue 
c ----- вычисление гран. условия для потенциала       
      if(dlev.ge.2) write(65,*) 'direct FFT'

      k1=km/2 
      
      do 12 l=1,lm1 
         s2=(amp+amg)*fig1(l) 
         do k=1,km 
            a(k)=s2                ! No changes - axial symmetry  
         enddo 
         call fftc(a,nk,-1,km) 
         do k=1,k1+1 
c****************************************************************	 
            if((k.gt.(rank*parakm)).and.
     =	       (k.le.((rank+1)*parakm))) then
        	ph(2*im1,k-rank*parakm,l)=dreal(a(k)) 
            endif
	 enddo 
         do k=2,k1 
            if(((k1+k).gt.(rank*parakm)).and.
     =	       ((k1+k).le.((rank+1)*parakm))) then
        	ph(2*im1,k1+k-rank*parakm,l)=dimag(a(k)) 
            endif
	 enddo 
c****************************************************************	 

   12 continue 

      ! EXPERIMENTAL TO SET NEUMANN BOUNDARY CONDITIONS 
      do 15 i=1,2*im1-1 
         s2=(amp+amg)*fig2(i) 
         do k=1,km 
            a(k)=s2      ! No changes due to initial symmetry
         enddo 
         call fftc(a,nk,-1,km) 
         do k=1,k1+1
            if((k.gt.(rank*parakm)).and.
     =	       (k.le.((rank+1)*parakm))) then
        	ph(i,k-rank*parakm,lm1)=dreal(a(k))
	    endif	 
         enddo 

         do k=2,k1
            if(((k+k1).gt.(rank*parakm)).and.
     =	       ((k+k1).le.((rank+1)*parakm))) then
                ph(i,k1+k-rank*parakm,lm1)=dimag(a(k)) 
	    endif
	 enddo 
   15 continue 
c**************************************************************
      do 115 i=1,2*im1-1 
         s2=(amp+amg)*fig2m1(i) 
         do k=1,km 
            a(k)=s2      ! No changes due to initial symmetry
         enddo 
         call fftc(a,nk,-1,km) 
         do k=1,k1+1
            if((k.gt.(rank*parakm)).and.
     =	       (k.le.((rank+1)*parakm))) then
        	ph(i,k-rank*parakm,lm)=dreal(a(k))
	    endif	 
         enddo 

         do k=2,k1
            if(((k+k1).gt.(rank*parakm)).and.
     =	       ((k+k1).le.((rank+1)*parakm))) then
                ph(i,k1+k-rank*parakm,lm)=dimag(a(k)) 
	    endif
	 enddo 
  115 continue 


      if(dlev.ge.2) write(65,*) 'FFT finished'
   
      do i=1,2*im1-1 
         fig2m1(i) = ph(i,1,lm) 

         do l=1,lm-1 
            ph(i,1,l)=ph(i,1,lm) 
            do k=2,parakm 
               ph(i,k,l)=0.0d0
            enddo 
         enddo 
      enddo

c**************************************************************

      s1=pi/km 
      do k=1,k1+1 
            si2seq(k)=(2.0d0*dsin((k-1.0d0)*s1)/hf)**2 
      enddo 
      
      do k=2,k1 
            si2seq(k+k1)=si2seq(k) 
      enddo 
c****************************************************      
      do k=1,km
         si2(k)=si2seq(k)
      enddo
c****************************************************      

c------------------------------------------------     
      call puas7(lit2)
c------------------------------------------------     
      
      omega=dsqrt(ro01*pi/4.0d0) 
      s1=4.0d0
      if(rank.eq.0) then 
          write(25,*) 'frequency =',omega 
          write(25,*) '   period =',s1 
      endif
      
      r0=rd1*dsqrt(2.0d0/3.0d0) 
      j=0 

      if(dlev.ge.1) write(65,*) 'distributing particles'

c*******************************************
      s1=pi*r0*r0/(jm*procs) 
      s2=sqrt(s1) 
      s3=2.0d0*s2*jm*procs/(r0*r0) 
c*******************************************
      s2=s2/2.0d0 
      s1=s1/pi 
      s8=0.0d0 
      s4=r0 
      sx=0.0d0 
      sy=0.0d0 

 2001 continue 
 
      s6=s3*(s4-s2) 
      k=0.5d0*(s6+s8+1.0d0) 
      k=2*k 
      s8=s8+s6-k 
      s7=abs(s4*s4-k*s1) 
      s7=sqrt(s7) 
      x=(s4+s7)/2.0d0 
      s4=s7 
      y=2.0d0*pi/k 
      xf0=y*g05cae(f9) 
      s9=rd1*dsqrt(1.0d0-(1.0d0-1.5d0*(x/rd1)**2)**(2.0d0/3.0d0))  
      
      s10=s9/hr+1.0d0 
      i=s10 
      s10=s10-i 
      er=(1.0d0-s10)*fr(i,2)+s10*fr(i+1,2) 
      us=-0.5d0*tau*er 
      vs=dsqrt(dabs(-er*s9-us*us)) 
      do 1009 j1=1,k 
      j=j+1 
c**********************************************
      if((j.gt.(jm*rank)).and.(j.le.(jm*(rank+1)))) then
          xr(j-jm*rank)=s9 
          xf(j-jm*rank)=y*(j1-1.0d0)+xf0
          vr(j-jm*rank)=g05dde(us,vd) 
          vf(j-jm*rank)=vs
	  qr(j-jm*rank)=qr1
          cs=cos(xf(j-jm*rank)) 
          sn=sin(xf(j-jm*rank)) 
          sx=sx+vr(j-jm*rank)*cs-vf(j-jm*rank)*sn 
          sy=sy+vr(j-jm*rank)*sn+vf(j-jm*rank)*cs 
      else
          t05dde=g05dde(us,vd)	  
      endif   
 1009 continue 
      if(j.lt.(jm*(rank+1))) goto 2001 

      call addvelocity(sx,sy)
      sx=sx/(jm*size) 
      sy=sy/(jm*size)

c**********************************************
      do 1010 j=1,jm 
      cs=cos(xf(j)) 
      sn=sin(xf(j)) 
      vx=vr(j)*cs-vf(j)*sn-sx 
      vy=vr(j)*sn+vf(j)*cs-sy 
      vr(j)=vx*cs+vy*sn 
      vf(j)=-vx*sn+vy*cs 
      qr(j)=qr1
 1010 continue 

c     tau1=0.5d0*tau 
c     do 11 i=2,im2 
c      do 11 k=1,km2 
c      s2=sigma(i,k) 
c      if((s2.lt.1d-8).or.(i.eq.im2)) then 
c         s3=0.0d0 
c      else 
c         s3=(pp(i+1,k)-pp(i-1,k))/(2.0d0*hr*s2) 
c      endif
c      s1=0.5d0*(fr(i,k)+fr(i-1,k))-s3 
c      ur(i,k)=-tau1*s1
c      
c      if((-rr(i)*s1-ur(i,k)**2).gt.0.0d0) then
c         uf(i,k)=dsqrt(-rr(i)*s1-ur(i,k)**2) 
c      else
c         uf(i,k) = 0.0d0	 
c      endif	 
      
   11 continue 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GAS
      do i = 1,im1+1
	   do k = 1,km2
	      pp(i,k) = 0
	      tt(i,k) = 0
	   enddo
	enddo


      omega=dsqrt(sigma0*pi/4.d0)
      tau1=0.5*tau
      do i=2,im2
         s2=sigma(i,3)
         if((s2.lt.1.e-8).or.(i.eq.im2)) then
            s3=0.d0
         else
            s3=0.9d0*omega*rr(i)
         endif
         do k=1,km2
            ur(i,k)=0.d0
            uf(i,k)=s3
         enddo
      enddo

      if(dlev.ge.2) write(65,*) 'preparing gas' 
      do 7 k=1,km2 
      ur(1,k)=0.0d0 
      uf(1,k)=-uf(2,k) 
    7 continue 

      i1=rd2/hr+1.d0
      do i=i1+1,im2
         pp(i,3)=0.d0
         tt(i,3)=0.d0
      enddo
      s1=-0.5d0*hr*sigma(i1,3)*
     =   (uf(i1,3)**2/rr(i1)+0.5d0*(fr(i1,3)+fr(i1-1,3)))
      tt(i1,3)=s1
      tt(i1-1,3)=3.d0*s1
c      call pr3(uf,'0. uf')
c      call pr3(sigma,'0. sigma')
c      call pr3(fr,'0. fr')
      do i=i1-1,2,-1
         tt(i-1,3)=tt(i+1,3)-2.d0*hr*sigma(i,3)*
     =   (uf(i,3)**2/rr(i)+0.5d0*(fr(i,3)+fr(i-1,3)))
      enddo
      do i=3,i1-2
         pp(i,3)=tt(i,3)+(tt(i+1,3)-2.d0*tt(i,3)+tt(i-1,3))/4.d0-
     =        (tt(i+2,3)-2.d0*tt(i,3)+tt(i-2,3))/16.d0
      enddo
      pp(2,3)=(tt(3,3)+2.d0*tt(2,3)+tt(1,3))/4.d0
c      pp(2,3)=tt(2,3)
      pp(1,3)=pp(2,3)
      pp(i1,3)=tt(i1,3)
      pp(i1-1,3)=tt(i1-1,3)
      do i=1,im1
         s1=pp(i,3)
         do k=1,km2
            pp(i,k)=s1
            if(sigma(i,k).gt.1.d-8) then
               tt(i,k)=s1/sigma(i,k)
            else
               tt(i,k)=0.d0
            endif
         enddo
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END GAS


c*********************************************************  
      jm=0.9*n_p/procs
c*********************************************************  
 
      call reduce(0,0d0,0) 

      if(dlev.ge.1) write(65,*) 'start procedure ends'

   14 continue 
      return 
      end 
c---------------------------------------------------------       
      subroutine final(m7,tx) 
      implicit none
      integer puasfile,n_p,procs,parakm,im3
      integer dlev
      integer rank,size 
      real*8  cx,cy,dx,dy

      integer zh,fsave,fmom,fcen,fgas,fen,ffall,pot3d,fwin,
     =        pard,parfour,parv,gasd,gasv,rotc,rdisp,pdisp,dynflag,
     =        rvel,rflow
      
      include 'pg2.par' 
      include 'para_pg2.par' 
      include 'auxillary.par'

c************************************
      integer m7,evn
c************************************
      	
      real*8 fi(im1+1,km2),rr(2*im1),si2(km),f(100),tx 
      real*8 ph(2*im1,2*parakm,lm1),sigma(im1+1,km2) 
	  real*8 ur(im1+1,km2),uf(im1+1,km2),pp(im1+1,km2) 
      real*8 xr(n_p),xf(n_p),vr(n_p),vf(n_p),qr(n_p)
      integer*4 ie(km),lep(2*parakm) 
      real*8 qr1,ams,aen,aim,tau,hr,hf,hz,threshold
      integer im2,nt,n7,jm,k,i,l,j,ml,n8,icen,nb
      integer prev,next,crank,harprocs,group(procs),harmonics(procs)
      integer lodge(procs,2*parakm),workload(procs),layers(km),up(procs)
      integer down(procs)
      character*14 tr,tr1 
      common/a/tau,hr,hf,hz,rr 
      common/a1/f,ml,n7,n8 
      common/ntc/nt
      common/b/ur,uf,pp 
      common/d/sigma 
      common/e/fi,ams,aen,aim 
      common/e1/jm,xr,xf,vr,vf
      common/e2/qr1,qr
      common/g1/ph 
      common/i/icen 
      common/si/si2 
      common/ts/ie,lep
c*************************************      
      common/lod/lodge
      common/wor/workload
      common/thd/threshold
      common/lay/layers
      common/har/harmonics
      common/harbounds/up,down
      common/grp/prev,next,crank,harprocs,group

      common/window/nb
      common/para1/rank,size
       
c      call CrytBegin
      
      evn = m7 - (m7/2)*2
      write(tr1,'(i2.2,i1.1)') rank,evn 
      tr='ptot'//tr1 
      tr=tr(1:7)//'.dat'       
c*************************************      

      open(10,file=tr,form='unformatted') 
      im2=im1+1 
      f(1)=tx 
      f(3)=nt 
      f(12)=qr1 
      f(16)=n7 
      f(22)=ams 
      f(23)=aen 
      f(24)=aim 
      f(25)=jm

      f(31) = cx
      f(32) = cy
      f(33) = dx
      f(34) = dy

      f(35) = threshold
      f(36) = prev
      f(37) = next
      f(38) = crank
      f(39) = harprocs

      write(10) (f(i),i=1,100) 

c     BALANCING DATA
      write(10) (harmonics(k),k=1,procs)
      write(10) ((lodge(i,k),i=1,procs),k=1,2*parakm)
      write(10) (workload(k),k=1,procs)
      write(10) (layers(k),k=1,km) 
      write(10) (up(k),k=1,procs)
      write(10) (down(k),k=1,procs)
      write(10) (group(k),k=1,procs)
       
      write(10) (ie(k),k=1,km) 
      write(10) (lep(k),k=1,2*parakm) 
      write(10) (rr(i),i=1,2*im1) 
      write(10) (si2(k),k=1,km) 
      write(10) (((ph(i,k,l),i=1,2*im1),k=1,2*parakm),l=1,lm1) 
      write(10) ((sigma(i,k),i=1,im2),k=1,km2) 
      write(10) ((pp(i,k),i=1,im2),k=1,km2) 
      write(10) ((ur(i,k),i=1,im2),k=1,km2) 
      write(10) ((uf(i,k),i=1,im2),k=1,km2) 
      write(10) (xr(j),j=1,jm) 
      write(10) (xf(j),j=1,jm) 
      write(10) (vr(j),j=1,jm) 
      write(10) (vf(j),j=1,jm)
      close(10)

c      call CrytEnd
       
      call synchronize

      if(rank.eq.0) then      
        open(89,file='check.txt',form='formatted') 
        write(89,1089) m7
 1089 	format(i4.4)      
        close(89)
      endif
       
    1 continue   
      return 
      end 

c	real*8 function amicro()
c	amicro = 0d0
c	end
c--------------------------------------------------------- 
      subroutine particle(m1) 
c      implicit real*8(a-h,o-z) 
      implicit none
      integer puasfile,n_p,procs,parakm,im3
      integer rank,size,dlev

c      include 'mpif.h' 
      include 'pg2.par' 
      include 'para_pg2.par' 

      real*8 f(100) 
      real*8 ro(im1+1,km2),vvr(im1+1,km2),vvf(im1+1,km2)  
      real*8 ur(im1+1,km2),uf(im1+1,km2),sigma(im1+1,km2),pp(im1+1,km2) 
      real*8 fr(im1+1,km2),ff(im1+1,km2),rr(2*im1) 
      real*8 xr(n_p),xf(n_p),vr(n_p),vf(n_p),qr(n_p)
      real*8 zm,pi,eps,rcen,hr,x,vrj,vfj,tau1,tau,s1,s2,s4,er,ef,
     =s5,us,tr,vs,x1,y1,s,sn,cs,u1,v1,af,y,t,hf,rd1,s3,rm,qr1,
     =hz,rd2
      real*8 dabs,dsqrt,ur1,uf1,qr2
      integer im2,im,km1,lm,icen,k,i,j,jm,k1,i1,m1,nt,ml,n7,n8,lit2
      
      character*4 st
      common/a/tau,hr,hf,hz,rr 
      common/a1/f,ml,n7,n8 
      common/ntc/nt
      common/b/ur,uf,pp 
      common/d/sigma 
      common/d1/ro 
      common/d2/tr,vvr,vvf 
      common/e1/jm,xr,xf,vr,vf
      common/e2/qr1,qr
      common/g/fr,ff 
      common/i/icen
      common/a2/zm,rd2,lit2,st 
      common/crm/rm,rd1
      common/para1/rank,size
      
      real*8 ens3,ens11,ens12,ens4,ens17,tjm,pars3,pars4,pars8
      common/pen/ens3,ens11,ens12,ens4,ens17,tjm
      common/rep/pars3,pars4,pars8
      common/vmax/ur1,uf1
      
c****************************************************
      common/parapart/partTime
      real*8 partTime,interval
             
      if(dlev.ge.1) write(65,*) 'starting particle ',nt
      call settime
c****************************************************       
      pi=3.14159265358979d0 
      eps=1d-8 
      im2=im1+1 
      im=im1-1 
      km1=km2-1 
      lm=lm1-1 
c      rm=f(7) 
      rcen=(icen-1.0d0)*hr 

      do 1 i=1,im2 
      do 1 k=1,km2 
      ro(i,k)=0.0d0 
      vvr(i,k)=0.0d0 
      vvf(i,k)=0.0d0 
    1 continue 

      do 2 j=1,jm 

      x=xr(j) 
      y=xf(j)
      qr2=qr(j) 
      vrj=vr(j) 
      vfj=vf(j) 
      t=0.0d0 
   77 tau1=tau-t 
      s1=tau1*dabs(vfj)-0.5d0*hf*x 
      if(s1.gt.0.) then 
         s1=0.5d0*hf*x/abs(vfj) 
         tau1=max(s1,0.01d0*tau) 
      endif 
      t=t+tau1 
       
      s1=x/hr+1.5d0 
      i=s1 
      s1=s1-i 
      s3=x/hr+1.0d0 
      i1=s3 
      s3=s3-i1 
      s2=y/hf+1.5d0 
      k=s2 
      s2=s2-k 
      s4=y/hf+1.0d0 
      k1=s4 
      s4=s4-k1 
      
      er=(1.0d0-s3)*((1.0d0-s2)*fr(i1,k)+s2*fr(i1,k+1))+ 
     =s3*((1.0d0-s2)*fr(i1+1,k)+s2*fr(i1+1,k+1)) 
      ef=(1.0d0-s1)*((1.0d0-s4)*ff(i,k1)+s4*ff(i,k1+1))+ 
     =s1*((1.0d0-s4)*ff(i+1,k1)+s4*ff(i+1,k1+1)) 
      s3=(1.0d0-s1)*((1.0d0-s2)*ur(i,k)+s2*ur(i,k+1))+ 
     =s1*((1.0d0-s2)*ur(i+1,k)+s2*ur(i+1,k+1)) 
      s4=(1.0d0-s1)*((1.0d0-s2)*uf(i,k)+s2*uf(i,k+1))+ 
     =s1*((1.0d0-s2)*uf(i+1,k)+s2*uf(i+1,k+1)) 
      s5=(1.0d0-s1)*((1.0d0-s2)*sigma(i,k)+s2*sigma(i,k+1))+ 
     =s1*((1.0d0-s2)*sigma(i+1,k)+s2*sigma(i+1,k+1)) 
c      if(dlev.ge.1) write(65,*) '2 ',jm
      us=vrj+tau1*(er+tr*s5*(s3-vrj)) 
      vs=vfj+tau1*(ef+tr*s5*(s4-vfj)) 
      x1=x+tau1*us 
      y1=tau1*vs 
      s=dsqrt(x1*x1+y1*y1) 
      sn=y1/s 
      cs=x1/s 
c      if(dlev.ge.1) write(65,*) '3 ',jm
      u1=cs*us+sn*vs 
      v1=cs*vs-sn*us 
      af=y+dasin(sn) 
      if(s.gt.rm) then 
         u1=-u1 
         s=2.0d0*rm-s 
      endif 
      if(af.lt.0.0d0) af=af+2d0*pi 
      if(af.gt.2.0d0*pi) af=af-2.0d0*pi 
c      if(dlev.ge.1) write(65,*) '4 ',jm
      x=s 
      y=af 
      vrj=u1 
      vfj=v1 
      if(abs(tau-t).gt.eps) goto 77 
      xr(j)=s 
      xf(j)=af 
      vr(j)=u1 
      vf(j)=v1
        
      if (ur1.lt.dabs(u1)) ur1=dabs(u1)
      if (uf1.lt.dabs(v1)) uf1=dabs(v1)
      if(s.lt.rcen) goto 2 
      s1=s/hr+1.5d0 
      i=s1 
      s1=s1-i 
      s2=af/hf+1.5d0 
      k=s2 
      s2=s2-k 
      s3=qr2*(1.0d0-s1)*(1.0d0-s2) 
      s4=s3*u1 
c      if(dlev.ge.1) write(65,*) '5 ',jm

 1075 format('pp ',i6,2e40.30)
      s5=s3*v1 
          ro(i,k)=ro(i,k)+s3 
          vvr(i,k)=vvr(i,k)+s4 
          vvf(i,k)=vvf(i,k)+s5 
      s3=qr2*(1.0d0-s1)*s2 
      s4=s3*u1 
      s5=s3*v1 
          ro(i,k+1)=ro(i,k+1)+s3 
          vvr(i,k+1)=vvr(i,k+1)+s4 
          vvf(i,k+1)=vvf(i,k+1)+s5 
      s3=qr2*s1*(1.0d0-s2) 
      s4=s3*u1 
      s5=s3*v1 
          ro(i+1,k)=ro(i+1,k)+s3 
          vvr(i+1,k)=vvr(i+1,k)+s4 
          vvf(i+1,k)=vvf(i+1,k)+s5 

      s3=qr2*s1*s2 
      s4=s3*u1 
      s5=s3*v1 
          ro(i+1,k+1)=ro(i+1,k+1)+s3 
          vvr(i+1,k+1)=vvr(i+1,k+1)+s4 
          vvf(i+1,k+1)=vvf(i+1,k+1)+s5 
	  
      
	 
c 1124 format('+ ',i3,e40.30)
c1125 format('+ ',e40.30)
c 1126 format('+ ',2e40.30)
c 1127 format('+ ',3e40.30)
            	  
    2 continue 

      do 3 i=1,im2 
      ro(i,2)=ro(i,2)+ro(i,km2) 
      ro(i,km1)=ro(i,km1)+ro(i,1) 
      ro(i,km2)=ro(i,2) 
      ro(i,1)=ro(i,km1) 
      vvr(i,2)=vvr(i,2)+vvr(i,km2) 
      vvr(i,km1)=vvr(i,km1)+vvr(i,1) 
      vvr(i,km2)=vvr(i,2) 
      vvr(i,1)=vvr(i,km1) 
      vvf(i,2)=vvf(i,2)+vvf(i,km2) 
      vvf(i,km1)=vvf(i,km1)+vvf(i,1) 
      vvf(i,km2)=vvf(i,2) 
      vvf(i,1)=vvf(i,km1) 
    3 continue 
      s1=(6.0d0*im+1.0d0)/(6.0d0*im-1.0d0) 
      s2=(6.0d0*icen-7.0d0)/(6.0d0*icen-5.0d0) 
      do 4 k=1,km2 
      ro(icen+1,k)=ro(icen+1,k)+s2*ro(icen,k) 
      ro(icen,k)=0.0d0 
      ro(im1,k)=ro(im1,k)+s1*ro(im2,k) 
      vvr(icen+1,k)=vvr(icen+1,k)+s2*vvr(icen,k) 
      vvr(icen,k)=0.0d0 
      vvr(im1,k)=vvr(im1,k)+s1*vvr(im2,k) 
      vvf(icen+1,k)=vvf(icen+1,k)+s2*vvf(icen,k) 
      vvf(icen,k)=0.0d0 
      vvf(im1,k)=vvf(im1,k)+s1*vvf(im2,k) 
    4 continue 

c*************************************************    
      call reduceparticles
      call particlenergy

      tjm = jm

      call addmatrix(ens3,ens11,ens12,ens4,ens17,
     =               tjm,pars3,pars4,pars8)

c*************************************************
      do k=1,km2 
         do i=2,im1 
            if(ro(i,k).gt.1.0d-8) then 
               s=ro(i,k) 
            else 
               s=1.0d0 
            endif   
            vvr(i,k)=vvr(i,k)/s 
            vvf(i,k)=vvf(i,k)/s 
            ro(i,k)=ro(i,k)/(hf*hr*rr(i)) 
	 enddo
      enddo	     

      do 6 k=1,km2 
      ro(1,k)=ro(2,k) 
      ro(im2,k)=ro(im1,k) 
    6 continue 

      call dens0(1,m1,ro) 

c****************************************************       
      partTime = interval()
c****************************************************       
      call adaptivemass(1)
      if(dlev.ge.1) write(65,*) 'finished particle'
      return 
      end 
      

      subroutine puas7(lit2) 
      implicit none 
c*************************************************
      integer zh,fsave,fmom,fcen,fgas,fen,ffall,pot3d,fwin,
     =        pard,parfour,parv,gasd,gasv,rotc,rdisp,pdisp,dynflag,
     =        rvel,rflow
      real*8 cx,cy,dx,dy

      integer puasfile,n_p,procs,parakm,im3,kmin,kmax
      integer rank,size 
c*************************************************
      	      
      integer im,im2,lm,i,k,l,km1,l3,lit,lit2 
      integer nt,k1,lk,gk,dlev
      include 'pg2.par' 
      include 'para_pg2.par' 
      include 'auxillary.par'

      real*8 fi(im1+1,km2),sigma(im1+1,km2), 
     =rr(2*im1),ro(im1+1,km2)
      
c***********************************************************     
      real*8 ph(2*im1,2*parakm,lm1),rh(im1,km),si2(km) 
      real*8 surph(im1+1,km),low(2*im1,lm1)
c***********************************************************     

      real*8 fr(im1+1,km2),ff(im1+1,km2) 
      real*8 pr(2*im1),amicro,pt1 
      real*8 s1,eps,hr,hf,hz,tau,omit,pi
      real*8 c1,c3,c4,ams0,hr2,hf2,eps1,ams,aim,aen,f(100)
      integer*4 ie(km),lep(2*parakm) 
      integer ml,n7,n8,getNegative
      complex*16 a(km) 
      common/a1/f,ml,n7,n8 
      common/a/tau,hr,hf,hz,rr 
      common/ntc/nt
      common/d/sigma 
      common/d1/ro 
      common/e/fi,ams,aen,aim 
      common/g/fr,ff 
      common/g1/ph 
      common/pf/pr,c1,c3 
      common/si/si2 
      common/ts/ie,lep 
      common/x2/hr2,hf2,ams0,eps,omit 
      common/rhdens/rh
c**********************************************************       
      common/para1/rank,size
      common/para2/surph 
      common/parapuas/puasTime
      real*8 puasTime
      integer lodge(procs,2*parakm)
      common/lod/lodge
      
      if(dlev.ge.1) write(65,*) 'starting Poisson ',nt
      print*,'starting Poisson ',nt
      
      
      pt1 = amicro()
c**********************************************************       
      eps1=eps/4.0d0 
      pi=3.14159265358979d0 
      lit=0 
      im2=im1+1 
      im=im1-1 
      km1=km2-1 
      lm=lm1-1 
      c1=1d0/hr**2 
      c3=1d0/hz**2 
      c4=4d0*pi/hz
      
c      call synchronize
      
c     0 means reading 0-th harmonic from file 
c     and skipping computations 
c     if not, working as usual 
      if ((ml.eq.0) .and. (rank.eq.0) .and. 
     =    (nt.eq.0) ) then

        open(66,file='low.hrm',form='unformatted') 

        read(66) ((low(i,l),i=1,2*im1),l=1,lm1) 

        close(66)

        do i = 1,2*im1
    	    do l = 1,lm1
		ph(i,1,l) = low(i,l)
	    enddo	
	enddo
	if(dlev.eq.2) write(65,*) 'read 1st harmonic ',ph(2,1,1)
	
      endif
      
c--------------------------- Фурье для плотности ---     
c****************************************
      
      k1 = km/2 
c****************************************

      do i=2,im1 
         do k=1,km 
            a(k)=ro(i,k+1)+sigma(i,k+1) 
         enddo 
         call fftc(a,nk,-1,km) 
         do k=1,k1+1 
            rh(i,k)=dreal(a(k)) 
         enddo 
         do k=2,k1 
            rh(i,k1+k)=dimag(a(k)) 
         enddo 
      enddo

c     if it is 0th step, of 1st started with given 0 harmonic, initializing balance arrays
c      if((nt.eq.0).or.((nt.eq.1).and.(ml.eq.0))) call setInitial
      
c------------------------------------- цикл по гармоникам      
      kmax = 2*parakm
      kmin = 1

      if ((rank.eq.0) .and. (nt.eq.0)) then

c     -1 means computation of only the lowest harmonics (k = 0,1)
          if (ml .eq. -1) then
    	    kmax = 1
	    kmin = 1
	  endif    

c     0 means reading 0th harmonic from a file, computing all the others
          if (ml .eq. 0) then
    	    kmax = 2*parakm
	    kmin = 2
	  endif    
      endif

      if(dlev.ge.2) write(65,*) 'Poisson prepared'
      print*,'Poisson prepared'
      

c     EVALUATING THE SPLIT HARMONIC (IF THERE is ONE)
c      lk = getNegative()

      do 12 lk=1,parakm 
	   gk = parakm*rank + lk
          print*,' axial harmonic ',lk,gk
         call poisson2D(gk,lk)
         stop
   12 continue                ! axial loop ends here

      if(dlev.ge.1) write(65,*) 'Axial loop finished'

c -------------------------- конец цикла по гармоникам      
      i=lit/(lm*km)+1 
      lit2=lit2+lit 
      l3=lit2/(lm*km)+1 
      
      if(rank.eq.0) then
        write(25,105) nt,i,l3 
      endif
      if(dlev.ge.5) write(65,*) 'A'
      	
  105 format('nt=',i5,'   number of iteration =',i4,'   total =',i6)    
c       
c***********************************************
c     saving the lowest harmonics
c     and exiting
      
      if (ml .eq. -1) then

        if(rank.eq.0) then

	    do i = 1,2*im1
		do l = 1,lm1
		    low(i,l) = ph(i,1,l)
		enddo
	    enddo

    	    open(66,file='low.hrm',form='unformatted') 

	    write(66) ((low(i,l),i=1,2*im1),l=1,lm1) 

	    close(66)
        endif
	
	call synchronize
	
	call finalize
	
	stop
      endif
c      
 5433 continue
      if(dlev.ge.5) write(65,*) 'B'

      call synchronize
      if(dlev.ge.5) write(65,*) 'B1'

      !gathering the sliced harmonics
      call gatherFi

      !m-m-making balance
      if(dynflag.eq.1) call balance(nt)
c      ***********************************************    
      if(dlev.ge.5) write(65,*) 'C'
              
      do i=2,im2 
         do k=1,k1+1 
            a(k)=dcmplx(surph(i,k),0d0) 
         enddo 
         do k=2,k1 
            a(k)=a(k)+dcmplx(0d0,surph(i,k1+k)) 
         enddo 
         do k=2,k1 
            a(k1+k)=dconjg(a(k1+2-k)) 
         enddo

	 call fftc(a,nk,1,km) 

         do k=1,km
	     fi(i,k+1)=dreal(a(k)) 
         enddo 
      enddo 
      
      do i=2,im2 
         fi(i,1)=fi(i,km1) 
         fi(i,km2)=fi(i,2) 
      enddo 
      do k=1,km2 
         fi(1,k)=fi(2,k) 
      enddo 

      do i=2,im1   
         s1=(ams0+ams)/((i-1d0)*hr)**2 
         do k=1,km2 
            fr(i,k)=(fi(i,k)-fi(i+1,k))/hr-s1
	 enddo
      enddo	     
      
      do k=1,km2 
         fr(1,k)=fr(2,k) 
      enddo 
      do 20 i=1,im2   
      do 20 k=1,km1 
      ff(i,k)=(fi(i,k)-fi(i,k+1))/(hf*rr(i)) 
   20 continue  
   
      if(dlev.ge.2) write(65,*) 
     =    'surface potential computed; gradient: ',fr(2,2),ff(2,2)

c************************************************
      puasTime = amicro() - pt1
c************************************************
c      call ClearMsg
      if(dlev.ge.1) write(65,*) 'Poisson finished'
      
      return 
      end 
c ---------------------------------------------- 
      subroutine pr1(b,i1,i2,j1,j2) 
      implicit none
      integer puasfile,n_p,procs,parakm,im3
      integer dlev
      include 'pg2.par' 
      include 'para_pg2.par' 

      real*8 vsp(7) 
      real*8 b(im1+1,km2) 
      integer isp(7),lm,n1,i1,n2,i2,ii1,i,l,j1,j2,j
      lm=lm1-1 
      n1=i1 
      n2=i1+6 
      if(n2.gt.i2) n2=i2 
 2003 continue 
 9001 format (/) 
      ii1=0 
      do 1001 i=n1,n2 
      ii1=ii1+1 
      isp(ii1)=i 
 1001 continue 
 
c      if(rank.eq.0) write(25,9002) (isp(i),i=1,ii1) 
      
 9002 format(/,7(8x,i2)) 
      do 1005 l=j1,j2 
      j=j1+j2-l 
      ii1=0 
      do 1003 i=n1,n2 
      ii1=ii1+1 
      vsp(ii1)=b(i,j) 
 1003 continue 
 
c      if(rank.eq.0) write(25,9004) j,(vsp(i),i=1,ii1) 

c 9004 format(1x,i3,7e10.3) 
 9004 format(1x,i3,7f10.5) 
 1005 continue 
c      if(rank.eq.0) write(25,9001) 
      if(n2.ge.i2) go to 2004 
      n1=n1+7 
      n2=n2+7 
      if(n2.gt.i2) n2=i2 
      go to 2003 
 2004 continue 
      return 
      end 

      subroutine density 
      implicit none
      integer puasfile,n_p,procs,parakm,im3
      integer dlev
      include 'pg2.par' 
      include 'para_pg2.par' 

      integer im2,im,km1,lm,i,k,j,jm,nt,ml,n7,n8
      real*8 x,y,s1,s2,hr,hf,hz,qr1,tau,qr2
      real*8 f(100) 
      real*8 ro(im1+1,km2),rr(2*im1) 
      real*8 xr(n_p),xf(n_p),vr(n_p),vf(n_p),qr(n_p) 
      common/a/tau,hr,hf,hz,rr 
      common/a1/f,ml,n7,n8 
      common/ntc/nt
      common/d1/ro 
      common/e1/jm,xr,xf,vr,vf
      common/e2/qr1,qr

      im2=im1+1 
      im=im1-1 
      km1=km2-1 
      lm=lm1-1 
      do 1 i=1,im2 
      do 1 k=1,km2 
      ro(i,k)=0d0 
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
          ro(i,k)=ro(i,k)+qr2*(1d0-s1)*(1d0-s2) 
          ro(i,k+1)=ro(i,k+1)+qr2*(1d0-s1)*s2 
          ro(i+1,k)=ro(i+1,k)+qr2*s1*(1d0-s2) 
          ro(i+1,k+1)=ro(i+1,k+1)+qr2*s1*s2 
    2 continue 
      do 3 i=1,im2 
      ro(i,2)=ro(i,2)+ro(i,km2) 
      ro(i,km1)=ro(i,km1)+ro(i,1) 
      ro(i,km2)=ro(i,2) 
      ro(i,1)=ro(i,km1) 
    3 continue 
      s1=(6d0*im+1d0)/(6d0*im-1d0) 
      do 4 k=1,km2 
      ro(2,k)=ro(2,k)-ro(1,k) 
      ro(im1,k)=ro(im1,k)+s1*ro(im2,k) 
    4 continue   
      do 5 k=1,km2 
      do 5 i=2,im1 
      ro(i,k)=ro(i,k)/(hr*hf*rr(i)) 
    5 continue 
      do 6 k=1,km2 
      ro(1,k)=ro(2,k) 
      ro(im2,k)=ro(im1,k) 
    6 continue 
      return 
      end 
c--------------------------------------------------------- 
      real*8 function press(ro) 
      implicit none
      real*8 cg,gamma,ro
      common/x4/cg,gamma 
      press=cg*ro**gamma 
      return 
      end 
c--------------------------------------------------------- 
      real*8 function ellip(a) 
      implicit none
      real*8 pi2,eps,s1,a,r,s,s2
      integer k
      
      pi2=1.570796326795d0 
      eps=1.e-8 
      s1=1d0/(1d0-a) 
      s=1d0 
      r=1d0 
      k=0 
    1 k=k+1 
      r=r*a*((2d0*k-1d0)/(2d0*k))**2 
      s=s+r 
      s2=s1*r 
      if(r.gt.eps) goto 1
      ellip=pi2*s 
      return 
      end 
c--------------------------------------------------------- 
      subroutine particlenergy
      implicit none
      
c*********************************************************       
      integer puasfile,procs,parakm,im3,dlev
c*********************************************************       

      integer j,n_p,i,k,im,lm,jm,icen,km1
      include 'pg2.par' 
      include 'para_pg2.par' 

      real*8 qr1,x,y,us,vs,rcen,t1,t2,t3,hr,hf,hz,tau,s
      real*8 pi,ams,aen,aim,qr2
      real*8 fi(im1+1,km2),rr(2*im1) 
      real*8 xr(n_p),xf(n_p),vr(n_p),vf(n_p),qr(n_p) 

      real*8 ens3,ens11,ens12,ens4,ens17,tjm

      common/a/tau,hr,hf,hz,rr 
      common/e/fi,ams,aen,aim 
      common/e1/jm,xr,xf,vr,vf
      common/e2/qr1,qr
      common/i/icen 
      common/pen/ens3,ens11,ens12,ens4,ens17,tjm

      pi=3.14159265358979d0 
      im=im1-1 
      km1=km2-1 
      lm=lm1-1 
      rcen=(icen-1d0)*hr 

      ens3=0d0 
      ens4=0d0 
      ens11=0d0 
      ens12=0d0 
      ens17=0d0 

      do 1 j=1,jm 
      s=xr(j) 
      y=xf(j)
      qr2=qr(j) 
      t1=s/hr+1.5d0 
      i=t1 
      t1=t1-i 
      t2=y/hf+1.5d0 
      k=t2 
      t2=t2-k 
c      if(j.eq.491166) write(65,*) 'pe 2',i,s,y
      t3=(1d0-t1)*((1d0-t2)*fi(i,k)+t2*fi(i,k+1))+ 
     +t1*((1d0-t2)*fi(i+1,k)+t2*fi(i+1,k+1)) 
      ens17=ens17+qr2*t3 
c      if(j.eq.491166) write(65,*) 'pe 3'
      us=vr(j) 
      vs=vf(j) 
      ens3=ens3+qr2*s*vs 
      ens4=ens4+qr2*(us**2+vs**2)
      x=s*dcos(xf(j)) 
      y=s*dsin(xf(j)) 
      ens11=ens11+qr2*x 
      ens12=ens12+qr2*y 
    1 continue 
      ens4=ens4/2d0 
      ens17=ens17/2d0
      
      tjm = jm 

      end
            
c--------------------------------------------------------- 
      subroutine energy(tx) 
      implicit none
      
c*********************************************************       
      integer puasfile,procs,parakm,rank,size,im3,dlev
      integer zh,fsave,fmom,fcen,fgas,fen,ffall,pot3d,fwin,
     =        pard,parfour,parv,gasd,gasv,rotc,rdisp,pdisp,dynflag,
     =        rvel,rflow
      real*8 taim,taen,cx,cy,dx,dy
c*********************************************************       

      integer n_p,i,k,im,km1,lm,jm,icen
      include 'pg2.par' 
      include 'para_pg2.par' 
      include 'auxillary.par'
      
      real*8 qr1,x,y,hr,hf,hz,tau,amp,amg, 
     =s,s1,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,tx,s15,s16, 
     =cg,gamma,c4,pi,s17,rcen,ams,aen,aim 
      real*8 fi(im1+1,km2),rr(2*im1) 
	  real*8 ur(im1+1,km2),uf(im1+1,km2),pp(im1+1,km2) 
	  real*8 sigma(im1+1,km2) 
      real*8 xr(n_p),xf(n_p),vr(n_p),vf(n_p),kappa,qr(n_p)

      real*8     ens3,ens11,ens12,ens4,ens17,tjm

      common/a/tau,hr,hf,hz,rr 
      common/b/ur,uf,pp 
      common/d/sigma 
      common/e/fi,ams,aen,aim 
      common/e1/jm,xr,xf,vr,vf
      common/e2/qr1,qr
      common/i/icen 
      common/x3/amp,amg 
      common/x4/gamma,kappa
      common/para1/rank,size
      common/pen/ens3,ens11,ens12,ens4,ens17,tjm

      if(dlev.ge.1) write(65,*) 'starting energy'
      pi=3.14159265358979d0 
      im=im1-1 
      km1=km2-1 
      lm=lm1-1 
      rcen=(icen-1d0)*hr 

      c4=1d0/(gamma-1d0) 
c       
	  s5=0d0 
	  s6=0d0 
      s7=0d0 
      s8=0d0 
      s9=0d0 
      s15=0d0 
      s16=0d0 
      do 10 i=icen+1,im1
	  do 10 k=2,km1 
      s=hf*(k-1.5d0) 
      x=rr(i)*dcos(s) 
      y=rr(i)*dsin(s) 
      s=rr(i)*sigma(i,k) 
      s5=s5+s                                 ! massa gaza 
      s6=s6+s*(ur(i,k)**2+uf(i,k)**2)         ! kin.en. gaza 
      s7=s7+s*rr(i)*uf(i,k)                   ! mom. gaza 
      s8=s8+s*x                               ! c.t. gaza 
      s9=s9+s*y                               ! c.t. gaza 
      s15=s15+rr(i)*pp(i,k)                   ! vnutr. en. gaza 
      s16=s16+s*fi(i,k)                       ! en. polya 
   10 continue 
	  s5=s5*hr*hf 
	  s6=0.5d0*s6*hr*hf 
      s7=s7*hr*hf 
      s8=s8*hr*hf 
      s9=s9*hr*hf 
      s15=c4*s15*hr*hf 
      s16=s16*hr*hf/2d0 
      
c*********************************************
c  gathering disributed values of invariants     
c  only those related to particles !!!!

      !ONLY AT INITIAL MOMENT: NORMALLY EVALUATED BY particle 
      if(tx.eq.0) call particlenergy
      
      s11 = ens11
      s12 = ens12
      s3  = ens3
      s4  = ens4
      s17 = ens17
      
      s11 = s11*jm        ! mass centre x
      s12 = s12*jm        ! mass centre y
      tjm  = jm           ! [M J+number of particles

      if(tx.eq.0) call gatherInv(s3,s11,s12,s4,s17,tjm,taim,taen)
      s11   = s11/tjm
      s12   = s12/tjm

c*********************************************
      
      s10=s3+s7+aim 
      s1=s4+s6+s15+s16+s17+aen  

      s13=(s8+s11)/(amg+amp) 
      s14=(s9+s12)/(amg+amp) 
      if(amg.eq.0d0) then 
         s8=0d0 
         s9=0d0 
      else 
         s8=s8/amg 
         s9=s9/amg 
      endif 
      if(amp.eq.0.) then 
         s11=0d0 
         s12=0d0 
      else 
         s11=s11/amp 
         s12=s12/amp 
      endif

      
      if(rank.eq.0) then 
          write(25,101) s5 
          write(25,111) s8,s9 
          write(25,112) s11,s12 
          write(25,113) s13,s14 
          write(25,121) s7 
          write(25,122) s3 
          write(25,123) s10 
          write(25,131) s6 
          write(25,132) s15 
          write(25,133) s4 
          write(25,134) s16 
          write(25,135) s17 
	  write(25,136) s1 
      endif
      	  
  101 format('     massa of gas =',e20.12)     
  111 format('     gas centre of gravity (x,y)=',2e14.6)     
  112 format('particle centre of gravity (x,y)=',2e14.6)     
  113 format('   total centre of gravity (x,y)=',2e14.6)     
  121 format('     moment of gas impulse =',e20.12)     
  122 format('moment of particle impulse =',e20.12)     
  123 format('   total moment of impulse =',e20.12)     
  131 format('     kinetic gas energy =',e20.12)     
  132 format('    inner energy of gas =',e20.12) 
  133 format('particle kinetic energy =',e20.12) 
  134 format('       gas field energy =',e20.12) 
  135 format('  particle field energy =',e20.12) 
  136 format('           total energy =',e20.12) 

      if(rank.eq.0) then 
          if(fmom.eq.1) write(31,155) tx,s5,s7,s3,s10 
	  if(fcen.eq.1) write(32,156) tx,s8,s9,s11,s12,s13,s14 
          if(fgas.eq.1) write(33,157) tx,s6,s15,s16 
    	  if(fen.eq.1)  write(34,158) tx,s4,s17,s1 
      endif
          	  
  155 format(f10.4,4e14.6) 
  156 format(f10.4,6e14.6) 
  157 format(f10.4,3e14.6) 
  158 format(f10.4,3e14.6) 
      if(dlev.ge.1) write(65,*) 'finished energy'
      return 
      end 
c--------------------------------------------------------- 
      subroutine dens0(n,m1,ro) 
      implicit none
      integer puasfile,n_p,procs,parakm,im3
      integer dlev
      include 'pg2.par' 
      include 'para_pg2.par' 
      
      integer lm,km1,k,k1,n,i,m1,m2,m
      real*8 s,hr,hf,g,a,b,tau,hz
      real*8 rr(2*im1),sn(km2),cs(km2),ro(im1+1,km2),pr(km2) 
      common/a/tau,hr,hf,hz,rr 
c      common/d1/ro
      
c 2596 format('@@ ',i3,e25.17)
c 2597 format('@ ',3e25.17)

c      write(65,2596) 0,ro(7,65)
      lm=lm1-1 
      km1=km2-1 
      k1=km/2 
      do k=1,km2 
         s=(k-1)*hf 
         sn(k)=dsin(s) 
         cs(k)=dcos(s) 
      enddo 
c      write(65,2596) 1,ro(7,65)
      g=0d0 
      a=0d0 
      b=0d0 
      do k=2,km1 
         g=g+ro(2,k) 
         a=a+ro(2,k)*cs(k) 
         b=b+ro(2,k)*sn(k) 
      enddo 
c      write(65,2596) 2,ro(7,65)
      g=g/km 
      a=2.*a/km 
      b=2.*b/km 
      if(n.eq.1) then 
         s=g*g-a*a-b*b 
         if(s.lt.0d0) then 
            s=g/dsqrt(a*a+b*b) 
            a=a*s 
            b=b*s 
         endif 
      else 
         g=0d0 
      endif 
c      write(65,2596) 3,ro(7,65)
      do k=1,km2 
         ro(2,k)=g+a*cs(k)+b*sn(k) 
      enddo 
c      write(65,2596) 4,ro(7,65)
      do k=1,km2 
         if(k.le.k1) then 
            ro(1,k)=ro(2,k+k1) 
         else 
            ro(1,k)=ro(2,k-k1) 
         endif 
      enddo 
c      write(65,2596) 5,ro(7,65)
c      write(65,2597) ro(7,1),0.0,ro(7,130)
      do i=3,im1 
         m2=m1*3.5d0/(2*i-3) 
         if(m2.lt.1) goto 1 
         do m=1,m2 
            do k=2,km1 
               pr(k)=0.25d0*(ro(i,k-1)+2d0*ro(i,k)+ro(i,k+1)) 
            enddo 
c	    if(i.eq.7) write(65,2597) pr(64),pr(65),pr(66)
            pr(1)=pr(km1) 
            pr(km2)=pr(2) 
            do k=2,km1 
c	       if((k.eq.65).and.(i.eq.7)) write(65,2596) -m,ro(7,65)
               ro(i,k)=0.25d0*(pr(k-1)+2d0*pr(k)+pr(k+1)) 
c	       if((k.eq.65).and.(i.eq.7)) write(65,2596) m,ro(7,65)
            enddo 
c	    if(i.eq.7) write(65,2597) ro(i,64),ro(i,65),ro(i,66)
            ro(i,1)=ro(i,km1) 
            ro(i,km2)=ro(i,2) 
         enddo 
      enddo 
c      write(65,2596) 6,ro(7,65)
    1 continue 
c      write(65,2596) 7,ro(7,65)
      return 
      end 
      
c--------------------------------------------------------- 
      real*8 function g05dde(a,b) 
      implicit real*8 (a-h,o-z) 
      real*8 a, b 
      real*8 store1, store2 
      real*8 half, one, t, u, v, w, x 
      integer n 
      real*8 d(41) 
      real*8 g05cae 
      common /cag05b/ store1, store2 
      data one /1.0d0/, half /0.5d0/ 
      data d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), d(9), 
     * d(10), d(11), d(12), d(13), d(14) /0.0d0,0.674489750196082d0, 
     * 1.150349380376008d0,1.534120544352546d0,1.862731867421652d0, 
     * 2.153874694061456d0,2.417559016236505d0,2.660067468617460d0, 
     * 2.885634912426757d0,3.097269078198785d0,3.297193345691964d0, 
     * 3.487104104114431d0,3.668329285121323d0,3.841930685501911d0/ 
      data d(15), d(16), d(17), d(18), d(19), d(20), d(21), d(22), 
     * d(23), d(24), d(25), d(26), d(27) /4.008772594168585d0, 
     * 4.169569323349106d0,4.324919040826046d0,4.475328424654204d0, 
     * 4.621231001499247d0,4.763001034267814d0,4.900964207963193d0, 
     * 5.035405969463927d0,5.166578119728753d0,5.294704084854598d0, 
     * 5.419983174916868d0,5.542594057802940d0,5.662697617459439d0/ 
      data d(28), d(29), d(30), d(31), d(32), d(33), d(34), d(35), 
     * d(36), d(37), d(38), d(39), d(40) /5.780439324478935d0, 
     * 5.895951216739571d0,6.009353565530745d0,6.120756285971941d0, 
     * 6.230260137989044d0,6.337957754553790d0,6.443934526538564d0, 
     * 6.548269367831731d0,6.651035379893011d0,6.752300431407015d0, 
     * 6.852127665896068d0,6.950575947916750d0,7.047700256664409d0/ 
      data d(41) /7.143552034352190d0/ 
      u = store1 
      do 20 n=1,39 
         if (u.gt.half) go to 40 
         u = u + u 
   20 continue 
      n = 40 
   40 t = d(n) 
      u = g05cae(x) 
   60 w = (d(n+1)-t)*u 
      v = w*(w*half+t) 
   80 u = g05cae(x) 
      if (v.le.u) go to 100 
      v = g05cae(x) 
      if (u.gt.v) go to 80 
      u = (v-u)/(one-u) 
      go to 60 
  100 u = (u-v)/(one-v) 
      if (u.gt.half) go to 120 
      store1 = u + u 
      g05dde = a + b*(w+t) 
      return 
  120 store1 = u + u - one 
      g05dde = a - b*(w+t) 
      return 
      end 
c--------------------------------------------------------- 
      real*8 function g05cae(x) 
      implicit real*8 (a-h,o-z) 
      real*8 x 
      integer ix,iy,iz 
      common /cag05a/ ix,iy,iz 
      ix = 171*(ix-(ix/177)*177) - 2*(ix/177) 
      iy = 172*(iy-(iy/176)*176) - 2*(iy/176) 
      iz = 170*(iz-(iz/178)*178) - 2*(iz/178) 
      if (ix.lt.0) ix = ix+30269 
      if (iy.lt.0) iy = iy+30307 
      if (iz.lt.0) iz = iz+30323 
      x=1. 
      ax=ix 
      ay=iy 
      az=iz 
      ai = ax/30269.0d0+ay/30307.0d0+az/30323.0d0 
      ii=ai 
      g05cae = ai-ii 
      return 
      end 
c--------------------------------------------------------- 
      block data 
      integer ix,iy,iz 
      real*8 normal,gamma 
      common /cag05a/ ix,iy,iz 
      common /cag05b/ normal,gamma 
      data ix/1/,iy/255/,iz/25555/ 
      data normal/1.0d0/,gamma/-1.0d0/ 
      end 
c--------------------------------------------------------- 
      integer function ich(st) 
      integer i,k 
      character*2 st 
      character*1 a 
      a=st(1:1) 
      i=ichar(a)-48 
      a=st(2:2) 
      k=ichar(a)-48 
      ich=10*i+k 
      return 
      end 

c--------------------------------------------------------- 

      subroutine reduceparticles 
      implicit none 

c*********************************************************       
      integer puasfile,procs,parakm,im3,dlev
c*********************************************************       

      integer n_p,km1,icen,jm,j,j1 
      include 'pg2.par' 
      include 'para_pg2.par' 

      real*8 qr1,us,vs,s,hr,hf,hz,tau
      real*8 pars3,pars4,pars8,rcen
      real*8 rr(2*im1),qr2 
      real*8 xr(n_p),xf(n_p),vr(n_p),vf(n_p),qr(n_p)
            
      common/a/tau,hr,hf,hz,rr 
      common/e1/jm,xr,xf,vr,vf
      common/e2/qr1,qr
      common/rep/pars3,pars4,pars8
      common/i/icen 
            

      rcen=(icen-1d0)*hr    
      km1=km2-1 
c       
      pars3=0d0 
      pars4=0d0 
      pars8=0d0 
      do 1 j=1,jm 
         s=xr(j) 
         if(s.lt.rcen) then 
            pars8=pars8+qr(j) 
            us=vr(j) 
            vs=vf(j) 
            qr2=qr(j)
            pars3=pars3+s*vs 
            pars4=pars4+us**2+vs**2 
            xr(j)=-2d0 
         endif  
    1 continue 

      j1=0 
      do 2 j=1,jm 
         if(xr(j).gt.0.) then 
            j1=j1+1 
            xr(j1)=xr(j) 
            xf(j1)=xf(j) 
            vr(j1)=vr(j) 
            vf(j1)=vf(j)
            qr(j1)=qr(j)
         endif 
    2 continue 
      jm=j1 

      pars3=pars3*qr2 
      pars4=pars4*qr2/2d0 
      

      end

c---------------------------------------------------------
      subroutine reduce(nt,tx,ffall) 
      implicit none 

c*********************************************************       
      integer puasfile,procs,parakm,rank,size,nt,im3
      integer ffall,dlev
      real*8 t1,t2,t3,t4,t5,t6,t7,t8,tx
c*********************************************************       

      integer n_p,i,k,km1,icen,jm
      include 'pg2.par' 
      include 'para_pg2.par' 

      real*8 qr1,hr,hf,hz,tau,amp,amg, 
     =s,s5,s6,s7,rcen,ams,aen,aim 
      real*8 fi(im1+1,km2),rr(2*im1) 
      real*8 ur(im1+1,km2),uf(im1+1,km2),pp(im1+1,km2) 
      real*8 sigma(im1+1,km2),ro(im1+1,km2) 
      real*8 xr(n_p),xf(n_p),vr(n_p),vf(n_p),qr(n_p)
      real*8 pars3,pars4,pars8
      common/a/tau,hr,hf,hz,rr 
      common/b/ur,uf,pp 
      common/d/sigma 
      common/d1/ro 
      common/e/fi,ams,aen,aim 
      common/e1/jm,xr,xf,vr,vf
      common/e2/qr1,qr
      common/i/icen 
      common/x3/amp,amg 
      common/para1/rank,size 
      common/rep/pars3,pars4,pars8
      
      if(dlev.ge.1) write(65,*) 'starting reduce'

      km1=km2-1 
c       
      s5=0d0 
      s6=0d0 
      s7=0d0
       
      rcen=(icen-1d0)*hr                  ! радиус центральной части 
      do 10 i=2,icen 
	  do 10 k=1,km2 
      s=rr(i)*sigma(i,k) 
      s5=s5+s                            !  massa gaza 
      s6=s6+s*(ur(i,k)**2+uf(i,k)**2)    !  energy gaza 
      s7=s7+s*rr(i)*uf(i,k)              !  impul's gaza 
      sigma(i,k)=0d0                      !  vnutri obolochki 
      ro(i,k)=0d0 
      ur(i,k)=0d0 
c      uf(i,k)=0d0 
   10 continue 
	  s5=s5*hr*hf 
	  s6=0.5d0*s6*hr*hf 
      s7=s7*hr*hf
c       

c********************************************************      
      if (nt.eq.0) call reduceparticles
              
      t1 = pars8 ! mass
      t2 = pars4 ! kynetic energy 
      t3 = pars3 ! momentum

      if (nt.eq.0) call gatherInv(t1,t2,t3,t4,t5,t6,t7,t8)

      ams=ams+t1+s5               !  massa Солнца 
      aen=aen+t2+s6               !  energy Солнца 
      aim=aim+t3+s7               !  moment impulsa Солнца 
c********************************************************      

c     writing the quantities accumulated in the centre
 1155 format(f10.4,3e14.6) 
      if((ffall.eq.1).and.(rank.eq.0)) then
        write(35,1155) tx,ams,aen,aim
      endif	

      
  100 format('n=',i4,'  jm=',i6,'  ams,aen,aim=',3e10.3) 
  
      if(dlev.ge.1) write(65,*) 'finished reduce'

      return 
      end 
     
c++++++++++++++++++++++++++++++++++++++++++++++
      subroutine fftc(a,n,isi,np) 
      implicit none 
      integer np,isi,i,nn,j,m,mm,ii,n,isign
      real*8  pi2,th
      complex*16 a(np),t,w,w1 
      
      do 7 i=1,np 
      if(isi.gt.0) a(i)=a(i)/np 
    7 continue 
      pi2=8.d0*datan(1.d0) 
      nn=np 
      j=1 
      do 3 i=1,nn 
      if(i.ge.j) go to 1 
      t=a(j) 
      a(j)=a(i) 
      a(i)=t 
    1 m=nn/2 
    2 if(j.le.m) go to 3 
      j=j-m 
      m=m/2 
      if(m.ge.1) go to 2 
    3 j=j+m 
      mm=1 
    4 if(mm.ge.nn) return 
      ii=2*mm 
c      th= pi2/dsign(ii,n*isi) 
      th= pi2/isign(ii,n*isi) 
      w1=dcmplx(-2.0d0*dsin(th/2)**2,dsin(th)) 
      w=1 
      do 6 m=1,mm 
      do 5 i=m,nn,ii 
      t=w*a(i+mm) 
      a(i+mm)=a(i)-t 
      a(i)=a(i)+t 
    5 continue 
      w=w1*w+w 
    6 continue 
      mm=ii 
      go to 4 
      end 

	real*8 function amicro()
	amicro = 0d0
	end


 
