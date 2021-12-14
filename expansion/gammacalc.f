      !
      !  Purpose: 
      !    This subroutine calculate the gamma function numerically, as it is needed for
      !    calculation of the particle diameters
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !


      subroutine gammacalc(visindex,val)
      use param
      implicit none
      
      ! defining parameters
      double precision, intent(in) :: visindex
      ! result
      double precision, intent(out) :: val
      
      ! exponent, sum, point of evaluation, dx for area calculation
      double precision :: ex, su, x, dx
      
      ! number of integration steps, loop index
      integer :: n, l
      
      
      ! calculate exponent for integration
      ex=5.d0/2.d0-visindex
      
      ! define number of integration steps
      n=200000000
      
      ! define length of integration steps
      dx=ex/dble(n)*24.d0
      
      ! integral sum
      su=0.d0
      ! loop over integration steps
      do l=1,n
        ! calculate position
        x=dble(l)*dx
        ! evaluate function at x and add to sum
        su=su+dexp((ex-1.d0)*dlog(x)-x)
      end do
      
      ! multiply by integration length to get the area
      val=su*dx
      
      
      
        
      
      
      
     
      end subroutine gammacalc
      

      
      
