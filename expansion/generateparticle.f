      !
      !  Purpose: 
      !    This subroutine generates the initial velocity of a particle
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !


      subroutine generateparticle(pt)
      use param
      implicit none
      
      ! particle type
      integer, intent(in) :: pt
      
      ! random numbers gaussian distributed
      real, dimension(1) :: rvec
      real, dimension(3) :: rvecwall
      double precision, dimension(3) :: vel
      integer :: indavp
      integer :: ix, iy, iz
      double precision :: phi, rad
      !integer :: pwall
      
      ! reset particle vector
      particle=0.d0
      
      ! generate a random point within the box
      call RANMAR(rvec, 1)
      rad=dsqrt(dble(rvec(1)))*radiusvalve
      call RANMAR(rvec, 1)
      phi=2*pi*dble(rvec(1))
      particle(1)=0.d0
      particle(2)=rad*dcos(phi)+disty/2.d0
      particle(3)=rad*dsin(phi)+distz/2.d0
      
      ! find the subcell
      ix=1
      iy=ceiling(particle(2)*ncelly/disty)
      iz=ceiling(particle(3)*ncellz/distz)
      
      ! save the subcell index
      particle(8)=dble(ix)
      particle(9)=dble(iy)
      particle(10)=dble(iz)
      


      
      call generatewallvel(pt,1,vel,1)


      
      ! save this velocity to as the particle velocity 
      particle(4:6)=vel
      ! save particle type
      particle(7)=dble(pt)
      
      ! set switch to 0, it should not record yet
      switchCounter=0
      parout=0
      
      ! generate random dt offset
      call RANMAR(rvec,1)
      dtprintoffset=dble(rvec(1))*dtprint(pt)
      

      
      end subroutine generateparticle
      

      
      
