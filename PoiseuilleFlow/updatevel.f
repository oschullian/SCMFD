      !
      !  Purpose: 
      !    This subroutine updates the avpvel whenever called
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  18/11/19   Otto Schullian    
      !


      subroutine updatevel
      use param
      implicit none
      
      ! particle type
      integer :: ptype
      ! subcell index in all directions
      integer :: ix,iy,iz
      

      ! retrieve particle type
      ptype=nint(particle(7))
      ix=nint(particle(8))
      iy=nint(particle(9))
      iz=nint(particle(10))
      
      ! calculate new position
      avtpos3(ptype,ix,iy,iz)=avtpos3(ptype,ix,iy,iz)+1
      
      ! check if it not bigger than the vector length
      if (avtpos3(ptype,ix,iy,iz)>nvelchange) 
     &   avtpos3(ptype,ix,iy,iz)=1
      
      ! add one to the number of particles to exchange at the next updating step,
      ! it is the minimum of the previous plus one, or the maximum number to exchange
      avtexchange3(ptype,ix,iy,iz)=
     &  min(avtexchange3(ptype,ix,iy,iz)+1,nvelchange)

      ! save velocity for the updating step
      avpvelinter3(avtpos3(ptype,ix,iy,iz),ptype,:,ix,iy,iz)
     & =particle(4:6)
      
      
      end subroutine updatevel
      

      
      
