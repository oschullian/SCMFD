      !
      !  Purpose: 
      !    This subroutine calculates the cummulative distribution function 
      !    of a gaussian probability distribution. This is needed to get 
      !    gaussian distributed 
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  17/02/18   Otto Schullian    
      !


      subroutine calccdf
      use param
      implicit none
      
      ! loop index
      integer :: k,l
      !integer :: k, l
      
      ! velocity
      double precision, dimension(ncdf) :: dummy_vcdf
      ! cdf of vel distribution function
      double precision, dimension(ncdf) :: dummy_cdf
      ! for the if-clause to find the corresponding dummy_cdf for a given cdf
      double precision :: testval
      ! distances to two points for linear interpolation
      double precision :: dx1, dx2 
      
      ! first we calculate a dummy cdf with equidistant velocities, 
      ! we will invert this as it is more convenient to have the cdf equidistant
      ! and the velocity at varying distances, the cdf is calculated in the unit sigma
      ! this way this one cdf can be used for all temperatures equally 
      
      ! initialise the arrays
      dummy_vcdf=0.d0
      dummy_cdf=0.d0
      
      ! first calculate the xvalues and f(x)
      do k=1,ncdf
        ! the velocity is linear with k and for k=1 it gives -8 and for k=ncdf
        ! it is 8 (8 * sigma is simply large enough for the erf
        ! to be close to 0 or 1, respectively
        dummy_vcdf(k)=8.d0*dble(k-1)/dble(ncdf-1)
     &   -8.e0*dble(k-ncdf)/dble(-ncdf+1)
        ! the corresponding error function
        dummy_cdf(k)=0.5d0*(1.d0-derf(-dummy_vcdf(k)/dsqrt(2.d0)))
      end do
      
      ! set the boundaries explicitly to 0 and 1 
      dummy_cdf(1)=0.d0
      dummy_cdf(ncdf)=1.d0
  
      ! now we invert this by defining the cdf equidistantly between 0 and 1 and find
      ! v with a linear interpolation
      
      ! loop over all elements on cdf
      do k=1,ncdf
        ! calculate cdf
        cdf(k)=dble(k-1)/dble(ncdf-1)
        
        
        do l=1,ncdf-1
           ! find two points on the dummy_cdf, where the cdf(k) is in between (then
           ! the testval will be negative
           testval=(dummy_cdf(l)-cdf(k))*(dummy_cdf(l+1)-cdf(k))
        
           if (testval.le.0.d0) then
            
             ! calculate the explicit distance of cdf(k) to both points 
             dx1=(cdf(k)-dummy_cdf(l))
             dx2=(dummy_cdf(l+1)-cdf(k))
             
             ! make sure that it is not zero, otherwise two different v have the same 
             ! dummy_cdf and we can not invert 
             if (dx1+dx2==0.d0) then
               write(*,*) 'cdf problem'
             end if
             
             ! calculate the corresponding velocity 
             vcdf(k)=dummy_vcdf(l)*dx2/(dx1+dx2)+
     &           dummy_vcdf(l+1)*dx1/(dx1+dx2)
             ! exit the inner loop and move to the next cdf point
             exit
             
           end if 
        end do 
      end do
      
      
   
      ! calculate the difference between two point of cdf, this will be used to 
      ! find the position of a random number on it. 
      dcdf=cdf(2)-cdf(1)
      
  
      end subroutine calccdf
      

      
      