	
	subroutine isnany(t,i)
	implicit none
	real*8 t
	integer i

	if(t.gt.0) then
	   i = 0
	else   
	   if(t.le.0) then
	      i = 0
	   else
	      i = 111
	   endif            
	endif
	
	end
