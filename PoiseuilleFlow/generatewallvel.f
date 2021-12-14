      !
      !  Purpose: 
      !    This subroutine updates the avpvel whenever called
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !


      subroutine generatewallvel(pt,walltype,vel,index)
      use param
      implicit none
      
      integer, intent(in) :: pt
      integer, intent(in) :: walltype
      integer, intent(in) :: index
      double precision, dimension(3), intent(out) :: vel
      
      real, dimension(1) :: rvec1dim
      
      double precision, dimension(3) :: velpar
      double precision :: vel1
      
      
      ! generate three velocities that are normal distributed
      call normrnd(velpar,3)
      ! multiply with sigma, because it is in units of sigma
      vel=velpar*dsqrt(Ak*tempwall/mass(pt))
      
      ! generate a random number and calculate a vx*exp(-vx.^2/2sigma) corresponding 
      ! velcoity 
      call RANMAR(rvec1dim, 1)
      vel1=dble(rvec1dim(1))
      vel1=dsqrt(-2.d0*dlog(1.d0-vel1))
      vel1=vel1*dsqrt(Ak*tempwall/mass(pt))
      
      vel(index)=vel1
      
      ! determine the sign depending on whether it was the left or right wall
      if (walltype==2) vel(index)=-vel(index)
      
      
      
      end subroutine generatewallvel
      

      
      