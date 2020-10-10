c-----initializing MPI library---------------------------------
      subroutine initialize(ml)

      integer procs,dlev,n_p,parakm,puasfile,im3

c      include 'mpif.h'
  
      include 'pg2.par' 
      include 'para_pg2.par' 

      integer rank, size, ierr, ml,proctype(procs),i
      character*12 tr,tr1
      common/para1/rank,size
      common/ptype/proctype

      rank = 0
	size = 1 
c      call MPI_Init(ierr)

c      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

c      call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

      
     
      end

c----------------------------------------------------
      subroutine initTime
      
      real*8 commTot,puasTot,gasTot,partTot,splitTime
      real*8 commTime,puasTime,gasTime,partTime
      integer rank,size
            
      common/timetot/commTot,puasTot,gasTot,partTot
      common/parapuas/puasTime
      common/paracomm/commTime
      common/split/splitTime
      common/paragas/gasTime
      common/parapart/partTime
      common/para1/rank,size

      commTot = 0.
      puasTot = 0.
      gasTot  = 0.
      partTot = 0.
      
      commTime = 0.
    
      if(rank.eq.0) then
          open(99,file='time.txt',form='formatted')
      endif
      
      end

c-----------------------------------------------      
c-----------------------------------------------      
      subroutine printTime(nt,flag)
      implicit none
      
      integer flag,nt,rank,size
      real*8 commTot,puasTot,gasTot,partTot,splitTime
      real*8 commTime,puasTime,gasTime,partTime

      common/timetot/commTot,puasTot,gasTot,partTot
      common/parapuas/puasTime
      common/paracomm/commTime
      common/paragas/gasTime
      common/parapart/partTime
      common/para1/rank,size
      common/split/splitTime

      if(rank.ne.0) return

 3142 format(i5,6f12.3)     

      if(flag.eq.0) then

        commTot = commTot + commTime
        partTot = partTot + partTime
        puasTot = puasTot + puasTime
        gasTot  = gasTot  + gasTime
      
        write(99,3142) nt,puasTime,splitTime,commTime,gasTime,partTime,
     =                    puasTime+gasTime+partTime 
      else
 3143   format(f20.3,4f15.3)     

        write(99,*) 'Total:'
        write(99,3143) ,puasTot,commTot,gasTot,partTot,
     =                  puasTot+gasTot+partTot 

        close(99)
	close(61)   !AWFUL
      endif

      commTime = 0.
      
      end      

                   
      subroutine settime
      real*8 first,amicro
      common/paratime2/first 

      first = amicro()
      end
      
      real*8 function interval()
      real*8 first,tt,amicro
      common/paratime2/first 

      tt = amicro()
      interval = tt - first
      
      end
                     
c--------closing MPI-----------------------------------------
      subroutine finalize
      integer parakm,procs,im3	
      integer n_p,puasfile,dlev,ierr,rank,size

c      include 'mpif.h'
  
      include 'pg2.par' 
      include 'para_pg2.par' 
      common/para1/rank,size
      
c      call MPI_Finalize(ierr)
      
      if(dlev.ge.1) write(65,*) ' finalised'
      if(dlev.gt.0) close(65)
      
      end

c------------------------------------------------------------    
      
c     similar to MPI_allgather, but sends to SELECTED procs
c     the data to send/result; size of the GLOBAL data       
      subroutine allgather(data1,sz)
      implicit none

      integer parakm,procs,im3	
      integer n_p,puasfile,dlev

c      include 'mpif.h'
  
      include 'pg2.par' 
      include 'para_pg2.par' 

      real*8 data1(4*km*im1),locdata(4*km*im1/procs)
      integer i,sz,proctype(procs),im2,im4
      integer rank,size,p,stat,ierr
      
      common/para1/rank,size
      common/ptype/proctype

      end

c----linkining potential harmonics together---------------------------------------------------------     
      subroutine gatherFi
      implicit none

      integer parakm,procs,im3	
      integer n_p,puasfile,dlev

c      include 'mpif.h'
  
      include 'pg2.par' 
      include 'para_pg2.par' 

      real*8 surph(im1 + 1,km)
      real*8 ph(2*im1,2*parakm,lm1)
      real*8 resfi((im1 + 4)*2*km),locfi((im1 + 4)*2*km)
      integer rank, size, k, i, ierr, im2      
      real*8 commTime,tt,amicro
      integer layers(km),lodge(procs,2*parakm),ti,proc,lock,im4,im5
      integer*4 ie(km),lep(2*parakm)
      integer prev,next,crank,harprocs,group(procs)

      common/ts/ie,lep 
      common/para1/rank,size	
      common/para2/surph
      common/g1/ph
      common/paracomm/commTime
      common/lay/layers
      common/lod/lodge
      common/grp/prev,next,crank,harprocs,group

      if(dlev.ge.1) write(65,*) 'preparing slice data'
	do i = 1,im1+1
	   do k = 1,km
	      surph(i,k) = ph(i,k,1)
	   enddo
	enddo

      end

c----gathering invariant values---------------------------------------------------------     
      subroutine gatherInv(f1,f2,f3,f4,f5,f6,f7,f8)
      
      implicit none

c      include 'mpif.h'
  

      integer rank, size, ierr
      real*8 commTime
      real*8 f1,f2,f3,f4,f5,f6,invec(20),outvec(20)
      real*8 f7,f8,tt,amicro

      common/para1/rank,size	
      common/paracomm/commTime
      
 
      end

      subroutine synchronize
c      include 'mpif.h'
  
      integer ierr

c      call MPI_Barrier(MPI_COMM_WORLD,ierr) 	
      end
            
      subroutine broadcast(m7)
c      include 'mpif.h'
  
      integer ierr,m7
      real*8 t
      
      
      end
                  
c----linkining scattered density matrix---------------------------------------------------------     
      subroutine addmatrix(f1,f2,f3,f4,f5,f6,f7,f8,f9)	

      implicit none	
      integer procs,im3	
      integer parakm,puasfile,n_p,dlev
      
c      include 'mpif.h'
  
      include 'pg2.par' 
      include 'para_pg2.par' 

      real*8 ro(im1 + 1,km2),vvr(im1+1,km2),vvf(im1+1,km2)
      real*8 rol((im1 + 1)*km2*3 + 9),res((im1 + 1)*km2*3 + 9)
      real*8 tr
      real*8 f1,f2,f3,f4,f5,f6,f7,f8,f9
      real*8 commTime,tt,amicro
      integer rank, size, l, k, ierr

      common/paracomm/commTime
      common/para1/rank,size	
      common/d1/ro
      common/d2/tr,vvr,vvf 

      end

      subroutine addvelocity(x,y)	

c      include 'mpif.h'
  

      real*8 x,y,ts(2),tr(2)
      integer rank, size,ierr

      common/para1/rank,size
      
      end

c----linkining scattered cartesian density matrix-------
      subroutine addcart	

      integer procs,dlev,puasfile,parakm,n_p,im3	

c      include 'mpif.h'
  
      include 'pg2.par' 
      include 'para_pg2.par' 

      real*8 rocart(im3,im3),velx(im3,im3),vely(im3,im3)
      real*8 jmfull(im1),vres(im1),rol(im3*im3),res(im3*im3)
      integer rank, size, l, k, ierr      

      common/para1/rank,size	
      common/para3/rocart,velx,vely,jmfull
            

      end

c----linkining scattered cartesian density matrix-------
      subroutine addmeanvalue(vrad,vphi)	

      integer procs,dlev,parakm,puasfile,n_p,im3	

c      include 'mpif.h'
  
      include 'pg2.par' 
      include 'para_pg2.par' 

      real*8 vrad(im1),vphi(im1)
      real*8 rol(2*im1),res(2*im1)
      integer rank, size, l, ierr      

      common/para1/rank,size	
            
      end
            
