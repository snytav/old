c ----------------------------------------------
      subroutine gaz(ams)
      implicit real*8(a-h,o-z)
	  include 'pg2.par'
      real*8 sigma(im1+1,km2),ur(im1+1,km2),uf(im1+1,km2),
     =pp(im1+1,km2),fr(im1+1,km2),ff(im1+1,km2),
     =ro(im1+1,km2),vvr(im1+1,km2),vvf(im1+1,km2),
     =ur1(im1+1,km2),uf1(im1+1,km2),rr(2*im1),
     =ee(im1+1,km2),tt(im1+1,km2),tt1(im1+1,km2)
      real*8 r4(im1),r2(im1),r3(im1)
      real*8 uu(im1+1,km2),vv(im1+1,km2),sigma1(im1+1,km2)
c      real*8 al(km2),be(km2),ga(km2),ps(km2),qs(km2)
      real*8 dr(km2),du(km2),dv(km2),de(km2),dp(km2),dt(km2)
      real*8 kappa
      common/a/tau,hr,hf,hz,rr
      common/b/ur,uf,pp
      common/c/ur1,uf1
      common/d/sigma
      common/d1/ro
      common/d2/trk,vvr,vvf
      common/i/icen
      common/g/fr,ff
      common/x4/gamma,kappa
	common/temp/tt
c      write(25,*) 'begin gaz'
c      write(25,*) 'tau=',tau
c      write(25,*) 'hr,hf,hz=',hr,hf,hz
c      write(25,*) 'trk=',trk
c      write(25,*) 'gamma,kappa=',gamma,kappa
c      print*,'begin gaz'
      pi=3.14159265358979d0
      im2=im1+1
      km1=km2-1
      km1=km+1
      lm=lm1-1
      c1=tau/hr
      c2=tau/hf
      c3=c1/2
      c4=c2/2
      tau1=0.5d0*tau
      eps=1.e-8
      do k=1,km2
c         dr(k)=1.5d0*sigma(3,k)-0.5d0*sigma(4,k)
         dr(k)=sigma(3,k)
         du(k)=ur(3,k)/3.d0
         dv(k)=uf(3,k)/3.d0
c         dp(k)=1.5d0*pp(3,k)-0.5d0*pp(4,k)
         dp(k)=pp(3,k)
      enddo
      do i=3,im2
         do k=1,km2
            if(sigma(i,k).gt.eps) then
               tt(i,k)=pp(i,k)/sigma(i,k)
            else
               tt(i,k)=0.d0
            endif
         enddo
      enddo
      do k=1,km2
         tt(2,k)=tt(3,k)
         tt(1,k)=tt(2,k)
      enddo
      do i=3,im1
         do k=2,km1
            if(ur(i,k).gt.0.d0) then
               s1=tt(i,k)-tt(i-1,k)
            else
               s1=tt(i+1,k)-tt(i,k)
            endif
            if(uf(i,k).gt.0.d0) then
               s2=tt(i,k)-tt(i,k-1)
            else
               s2=tt(i,k+1)-tt(i,k)
            endif
            s3=1.d0/(gamma-1.d0)+
     =         c3*(rr(i+1)*ur(i+1,k)-rr(i-1)*ur(i-1,k))/rr(i)+
     =         c4*(uf(i,k+1)-uf(i,k-1))/rr(i)
            s4=tau*trk*ro(i,k)*((vvr(i,k)-ur(i,k))**2+
     =         (vvf(i,k)-uf(i,k))**2)-
     =         c1*s1*ur(i,k)/(gamma-1.d0)-
     =         c2*s2*uf(i,k)/(rr(i)*(gamma-1.d0))+
     =         tt(i,k)/(gamma-1.d0)
            tt1(i,k)=s4/s3
            if(tt1(i,k).lt.0.d0) tt1(i,k)=0.d0
         enddo
      enddo
      do k=2,km1
         tt1(2,k)=tt1(3,k)
         tt1(1,k)=tt1(2,k)
         tt1(im2,k)=tt1(im1,k)
      enddo
      do i=1,im2
         tt1(i,1)=tt1(i,km1)
         tt1(i,km2)=tt1(i,2)
      enddo
c      write(25,*) '02.  tt'
c      call pr1(tt,40,46,82,98)
c      write(25,*) '03.  pp'
c      call pr1(pp,40,46,82,98)
c      
c      call pr3(ur,'1. ur')
c-------------------------- begin of 1 step -------------
      do k=1,km2
         sigma(2,k)=dr(k)
         sigma(1,k)=dr(k)
         ur(2,k)=du(k)
         ur(1,k)=-du(k)
         uf(2,k)=dv(k)
         uf(1,k)=-dv(k)
         pp(2,k)=dp(k)
         pp(1,k)=dp(k)
      enddo
      do i=1,im2
         do k=1,km2
            uu(i,k)=sigma(i,k)*ur(i,k)
            vv(i,k)=rr(i)*sigma(i,k)*uf(i,k)
            ee(i,k)=pp(i,k)/(gamma-1.d0)+
     =              sigma(i,k)*(ur(i,k)**2+uf(i,k)**2)/2.d0
         enddo
      enddo
c      call priv(uu,'uu01.txt')
c      write(25,*) '1. tau=',tau
c      write(25,*) '1. hr=',hr
c      write(25,*) '1. c3=',c3
c      write(25,*) '1. tau1=',tau1
c      write(25,*) '1. trk=',trk
c      write(25,*) '1. rr'
c      write(25,100) (rr(i),i=1,23)
c  100 format(6e12.4)
c      call pr3(uu,'1. uu')
c      call pr3(uf,'1. uf')
c      call pr3(pp,'1. pp')
c      call pr3(fr,'1. fr')
      do i=2,im1
         do k=2,km1
            s=sigma(i,k)
            if(s.gt.eps) then
c               if(uu(i,k).gt.0.d0) then
c                  s3=pp(i+1,k)-pp(i,k)
c               else
c                  s3=pp(i,k)-pp(i-1,k)
c               endif
c               if(vv(i,k).gt.0.d0) then
c                  s4=pp(i,k+1)-pp(i,k)
c               else
c                  s4=pp(i,k)-pp(i,k-1)
c               endif
c               if(k.eq.3) then
c                  r4(i)=tau*s*uf(i,k)**2/rr(i)
c                  r2(i)=c3*(pp(i+1,k)-pp(i-1,k))
c                  r3(i)=tau1*s*(fr(i,k)+fr(i-1,k))
c               endif
               s1=uu(i,k)+tau*s*uf(i,k)**2/rr(i)-
     =               c3*(pp(i+1,k)-pp(i-1,k))+
     =               tau1*s*(fr(i,k)+fr(i-1,k))
c     =              +tau*trk*s*ro(i,k)*(vvr(i,k)-ur(i,k))
               s2=vv(i,k)-tau*s*ur(i,k)*uf(i,k)-
     =               c4*(pp(i,k+1)-pp(i,k-1))+
     =               tau1*s*(ff(i,k)+ff(i,k-1))
c     =              +tau*trk*s*ro(i,k)*rr(i)*(vvf(i,k)-uf(i,k)) 
               ee(i,k)=ee(i,k)-c3*(pp(i+1,k)*ur(i+1,k)-
     =            pp(i-1,k)*ur(i-1,k))-tau*pp(i,k)*ur(i,k)/rr(i)-
     =            (c4/rr(i))*
     =            (pp(i,k+1)*uf(i,k+1)-pp(i,k-1)*uf(i,k-1))+
     =            s*(tau1*ur(i,k)*(fr(i,k)+fr(i-1,k))+
     =            tau1*uf(i,k)*(ff(i,k)+ff(i,k-1)))
c     =           +tau*trk*s*ro(i,k)*(vvr(i,k)*(vvr(i,k)-ur(i,k))+
c     =            vvf(i,k)*(vvf(i,k)-uf(i,k)))
               if(ee(i,k).lt.0.) ee(i,k)=0.
               uu(i,k)=s1
               vv(i,k)=s2
			endif
		 enddo
      enddo
c      call pr4(r4,'2. r4')
c      call pr4(r2,'2. r2')
c      call pr4(r3,'2. r3')
c      call pr3(uu,'2. uu')
c      open(27,file='r423.dat',form='formatted')
c      do i=1,10
c         write(27,100) i,r4(i),r2(i),r3(i)
c  100    format(i4,3e12.4)
c      enddo
c      close(27)
c      pause 11
c      call priv(uu,'uu02.txt')
c-------------------------- end of 1 step -------------
      do i=1,im2
         do k=1,km2
            if(sigma(i,k).le.eps) then
               sigma1(i,k)=0.
               ur1(i,k)=0.
			   uf1(i,k)=0.
               ee(i,k)=0.
               ur(i,k)=0.
               uf(i,k)=0.
            else
               sigma1(i,k)=sigma(i,k)*rr(i)
               ur1(i,k)=uu(i,k)*rr(i)
               uf1(i,k)=vv(i,k)*rr(i)
               ee(i,k)=ee(i,k)*rr(i)
               ur(i,k)=uu(i,k)/sigma(i,k)
               uf(i,k)=vv(i,k)/(sigma(i,k)*rr(i))
            endif
         enddo
      enddo
      do k=1,km2
         sigma1(2,k)=0.
         ur1(2,k)=0.
         uf1(2,k)=0.
         ee(2,k)=0.
         ur(2,k)=0.
         uf(2,k)=0.
         sigma1(1,k)=0.
         ur1(1,k)=0.
         uf1(1,k)=0.
         ee(1,k)=0.
         ur(1,k)=0.
         uf(1,k)=0.
      enddo
c      pause 12
      do i=2,im1
         ur(i,1)=ur(i,km1)
         ur(i,km2)=ur(i,2)
         uf(i,1)=uf(i,km1)
         uf(i,km2)=uf(i,2)
      enddo
      do k=1,km2
         ur(1,k)=-ur(2,k)
         ur(im2,k)=-ur(im1,k)
         uf(1,k)=-uf(2,k)
         uf(im2,k)=uf(im1,k)
      enddo
c      call pr3(ur,'3. ur')
c      write(25,*) '04.  sigma'
c      call pr1(sigma,40,46,82,98)
c      write(25,*) '05.  ur'
c      call pr1(ur,40,46,82,98)
c      write(25,*) '06.  uf'
c      call pr1(uf,40,46,82,98)
c      write(25,*) '07.  ee'
c      call pr1(ee,40,46,82,98)
c      write(25,*) '08.  pp'
c      call pr1(pp,40,46,82,98)
c---------------------- 1 begin of 2 step -----------------      
	  do i=1,im2
		 do k=1,km2
            sigma(i,k)=0.d0
            uu(i,k)=0.d0
            vv(i,k)=0.d0
            tt(i,k)=0.d0             ! new ee
         enddo
      enddo
c---------------------- 2 begin of 2 step -----------------      
      amu=0.49995*hf
      aru=0.49995*hr
c      print*,'end of part1'
c      pause 13
      do i=2,im1                                     ! (1)
         do k=2,km1                                  ! (2)
            r=hr*(i-1.5)      ! nach. koord. chastitsy
            fi=hf*(k-1.5)     ! nach. koord. chastitsy
c opredelenie nach. skorostej
            t=0.
            tau1=tau
   77       s1=r/hr
            i1=s1+1.5
            s1=s1+1.5-i1
            s2=fi/hf
            k1=s2+1.5
            s2=s2+1.5-k1
            s4=(1.-s1)*((1.-s2)*sigma1(i1,k1)+s2*sigma1(i1,k1+1))+
     =         s1*((1.-s2)*sigma1(i1+1,k1)+s2*sigma1(i1+1,k1+1))
            if(s4.gt.eps) then
               ux=((1.-s1)*((1.-s2)*ur(i1,k1)*sigma1(i1,k1)+
     =                  s2*ur(i1,k1+1)*sigma1(i1,k1+1))+
     =         s1*((1.-s2)*ur(i1+1,k1)*sigma1(i1+1,k1)+
     =                  s2*ur(i1+1,k1+1)*sigma1(i1+1,k1+1)))/s4
               uy=((1.-s1)*((1.-s2)*uf(i1,k1)*sigma1(i1,k1)+
     =                  s2*uf(i1,k1+1)*sigma1(i1,k1+1))+
     =         s1*((1.-s2)*uf(i1+1,k1)*sigma1(i1+1,k1)+
     =                  s2*uf(i1+1,k1+1)*sigma1(i1+1,k1+1)))/s4
            else
               ux=0.
               uy=0.
            endif
c            print*,'ux,uy=',ux,uy
   76       r1=r+tau1*ux
            y=tau1*uy
            sn=y/r
            alpha=dasin(sn)
            s3=dabs(alpha)
            if(s3.gt.amu) then  
               tau1=0.5*tau1
               goto 76
            endif
            s3=dabs(r1-r)
            if(s3.gt.aru) then  
               tau1=0.5*tau1
               goto 76
            endif
            r=r1
            if(r.gt.rr(im1)) r=rr(im1)
            fi=fi+alpha
            if(fi.lt.0.) fi=fi+2.*pi
            if(fi.gt.2.*pi) fi=fi-2.*pi
            t=t+tau1
c            s=tau1
            tau1=tau-t
c            print 300,i,k,s,tau1,t,tau
c  300       format('i,k,tau1,tau1,t,tau=',2i4,4e12.4)
c            pause 41
            if(tau1.gt.eps) goto 77
c-------------------------------------               
            s1=r/hr
            i1=s1+1.5
            s1=s1+1.5-i1
            s1=s1*(1.d0-(1.d0-s1)*hr/rr(i1))       !  ****
            s2=fi/hf
            k1=s2+1.5
            s2=s2+1.5-k1
            sigma(i1,k1)=sigma(i1,k1)+sigma1(i,k)*(1.-s1)*(1.-s2)
            sigma(i1,k1+1)=sigma(i1,k1+1)+sigma1(i,k)*(1.-s1)*s2
            sigma(i1+1,k1)=sigma(i1+1,k1)+sigma1(i,k)*s1*(1.-s2)
            sigma(i1+1,k1+1)=sigma(i1+1,k1+1)+sigma1(i,k)*s1*s2
            uu(i1,k1)=uu(i1,k1)+ur1(i,k)*(1.-s1)*(1.-s2)
            uu(i1,k1+1)=uu(i1,k1+1)+ur1(i,k)*(1.-s1)*s2
            uu(i1+1,k1)=uu(i1+1,k1)+ur1(i,k)*s1*(1.-s2)
            uu(i1+1,k1+1)=uu(i1+1,k1+1)+ur1(i,k)*s1*s2
            vv(i1,k1)=vv(i1,k1)+uf1(i,k)*(1.-s1)*(1.-s2)
            vv(i1,k1+1)=vv(i1,k1+1)+uf1(i,k)*(1.-s1)*s2
            vv(i1+1,k1)=vv(i1+1,k1)+uf1(i,k)*s1*(1.-s2)
            vv(i1+1,k1+1)=vv(i1+1,k1+1)+uf1(i,k)*s1*s2
            tt(i1,k1)=tt(i1,k1)+ee(i,k)*(1.-s1)*(1.-s2)
            tt(i1,k1+1)=tt(i1,k1+1)+ee(i,k)*(1.-s1)*s2
            tt(i1+1,k1)=tt(i1+1,k1)+ee(i,k)*s1*(1.-s2)
            tt(i1+1,k1+1)=tt(i1+1,k1+1)+ee(i,k)*s1*s2
         enddo
      enddo
c      pause 14
c----------------------gran.usloviya         
	  do i=1,im2
         sigma(i,2)=sigma(i,2)+sigma(i,km2)
         sigma(i,km1)=sigma(i,km1)+sigma(i,1)
         sigma(i,km2)=sigma(i,2)
         sigma(i,1)=sigma(i,km1)
         uu(i,2)=uu(i,2)+uu(i,km2)
         uu(i,km1)=uu(i,km1)+uu(i,1)
         uu(i,km2)=uu(i,2)
         uu(i,1)=uu(i,km1)
         vv(i,2)=vv(i,2)+vv(i,km2)
         vv(i,km1)=vv(i,km1)+vv(i,1)
         vv(i,km2)=vv(i,2)
         vv(i,1)=vv(i,km1)
         tt(i,2)=tt(i,2)+tt(i,km2)
         tt(i,km1)=tt(i,km1)+tt(i,1)
         tt(i,km2)=tt(i,2)
         tt(i,1)=tt(i,km1)
      enddo
	  do k=1,km2
         sigma(2,k)=sigma(2,k)+sigma(1,k)
         sigma(im1,k)=sigma(im1,k)+sigma(im2,k)
         sigma(1,k)=0.
         sigma(im2,k)=0.
         uu(2,k)=uu(2,k)+uu(1,k)
         uu(im1,k)=uu(im1,k)+uu(im2,k)
         uu(1,k)=0.
         uu(im2,k)=0.
         vv(2,k)=vv(2,k)+vv(1,k)
         vv(im1,k)=vv(im1,k)+vv(im2,k)
         vv(1,k)=0.
         vv(im2,k)=0.
         tt(2,k)=tt(2,k)+tt(1,k)
         tt(im1,k)=tt(im1,k)+tt(im2,k)
         tt(1,k)=0.
         tt(im2,k)=0.
      enddo
c      call pr3(uu,'4. uu')
      do i=2,im1
		 do k=1,km2
            if(sigma(i,k).le.eps) then
               ur(i,k)=0.
               uf(i,k)=0.
               sigma(i,k)=0.
               ee(i,k)=0.
            else
               ur(i,k)=uu(i,k)/sigma(i,k)
               uf(i,k)=vv(i,k)/(sigma(i,k)*rr(i))
               sigma(i,k)=sigma(i,k)/rr(i)
               ee(i,k)=tt(i,k)/rr(i)
            endif
         enddo
      enddo
c      pause 15
      do k=1,km2
         sigma(1,k)=sigma(2,k)
         sigma(im2,k)=sigma(im1,k)
         ur(1,k)=-ur(2,k)
         ur(im2,k)=ur(im1,k)
         uf(1,k)=-uf(2,k)
         uf(im2,k)=uf(im1,k)
         ee(1,k)=ee(2,k)
         ee(im2,k)=ee(im1,k)
      enddo
c      call pr3(ur,'5. ur')
      call dens0(1,2,sigma)
      call dens0(2,2,ur)
      call dens0(2,2,uf)
      call dens0(1,2,ee)
c      call pr3(ur,'6. ur')
c      pause 16
c-------- new --- new --- new --- new --- new --- new ---      
      do i=2,im1
         do k=1,km2
            if(sigma(i,k).gt.eps) then
               pp(i,k)=sigma(i,k)*tt1(i,k)
               s1=2.d0*(ee(i,k)-pp(i,k)/(gamma-1.d0))/sigma(i,k)
               s2=ur(i,k)**2+uf(i,k)**2
               if((s1.lt.0.d0).or.(s2.le.0.d0)) then
                  s1=0.d0
                  pp(i,k)=(gamma-1.d0)*ee(i,k)
               endif
               s=dsqrt(s1/s2)
			   ur(i,k)=s*ur(i,k)
			   uf(i,k)=s*uf(i,k)
            else
               pp(i,k)=0.d0
               ur(i,k)=0.d0
               uf(i,k)=0.d0
			endif	
         enddo
      enddo
      do k=1,km2
         pp(1,k)=pp(2,k)
      enddo
c      call pr3(ur,'7. ur')
c      pause 17
c-------- end new --- end new --- end new --- end new ---      
c      do i=2,im1
c         do k=1,km2
c            s1=sigma(i,k)*(ur(i,k)**2+uf(i,k)**2)/2.
c            pp(i,k)=(gamma-1.)*(ee(i,k)-s1)
c            if(pp(i,k).lt.0.) then
c               s=dsqrt(ee(i,k)/s1)
c               pp(i,k)=0.
c               ur(i,k)=s*ur(i,k)
c               uf(i,k)=s*uf(i,k)
c            endif   
c         enddo
c      enddo
cc      do k=1,km2
c         pp(icen,k)=pp(icen+1,k)
cc         pp(2,k)=2.d0*pp(3,k)-pp(4,k)
cc         pp(1,k)=pp(2,k)
cc         pp(im2,k)=pp(im1,k)
cc      enddo
c      write(25,*) '11.  sigma'
c      call pr1(sigma,40,46,82,98)
c      write(25,*) '12.  ur'
c      call pr1(ur,40,46,82,98)
c      write(25,*) '13.  uf'
c      call pr1(uf,40,46,82,98)
c      write(25,*) '14.  ee'
c      call pr1(ee,40,46,82,98)
c      write(25,*) '15.  pp'
c      call pr1(pp,40,46,82,98)
c      print*,'end gaz'
c      write(25,*) 'end gaz'
      return
      end
