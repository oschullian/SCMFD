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
      !integer :: pwall
      


      
      ! reset particle vector
      particle=0.d0
      
      ! generate a random point within the box
      call RANMAR(rvecwall, 3)
      particle(1)=dble(rvecwall(1))*distx
      particle(2)=dble(rvecwall(2))*disty
      particle(3)=dble(rvecwall(3))*distz
      
      ! find the subcell
      ix=ceiling(rvecwall(1)*real(ncellx))
      iy=ceiling(rvecwall(2)*real(ncelly))
      iz=ceiling(rvecwall(3)*real(ncellz))
      
      ! save the subcell index
      particle(8)=dble(ix)
      particle(9)=dble(iy)
      particle(10)=dble(iz)
      
      ! generate a random index between 0 and nopcopos 
      call RANMAR(rvec, 1)
      indavp=ceiling(real(nopcopos3(pt,ix,iy,iz))
     &   *rvec(1))
      
      
      

      ! save this velocity to as the particle velocity 
      particle(4:6)=avpvel3(indavp,pt,:,ix,iy, iz)
      ! save particle type
      particle(7)=dble(pt)
      
      ! set switch to 0, it should not record yet
      switchCounter=0
      
      ! generate random dt offset
      call RANMAR(rvec,1)
      dtprintoffset=dble(rvecwall(1))*dtprint(pt)
      

      
      end subroutine generateparticle
      

      
      
