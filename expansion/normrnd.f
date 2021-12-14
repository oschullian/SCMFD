      !
      !  Purpose: 
      !    Generate len normal distributed random numbers
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !


      subroutine normrnd(rvec,len)
      use param
      implicit none
      
      integer, intent(in) :: len
      double precision, intent(out), dimension(len) :: rvec
      real, dimension(len) :: rvecrvec
      integer :: k, l
      double precision :: dx1, dx2
      
      ! generate len uniformly distributed numbers
      call RANMAR(rvecrvec, LEN)
      
      ! make them double precision
      rvec=dble(rvecrvec)
      
      ! loop over all the randomly generated numbers
      do k=1,len
      
      ! calculate the index in cdf
      l= floor(rvec(k)/dcdf) +1

      ! calculate the distances of to the two point around 
      
      dx1=(rvec(k)-cdf(l))
      dx2=(cdf(l+1)-rvec(k))
             
      if (dx1+dx2==0.d0) then
        write(*,*) 'normrnd problem'
       end if
             
             
   
       ! linearly interpolate between these numbers
       rvec(k)=vcdf(l)*dx2/(dx1+dx2)+
     &           vcdf(l+1)*dx1/(dx1+dx2)
        
        
       
      end do
    
   
      
      end subroutine normrnd
      
      
      
      