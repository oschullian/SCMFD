      !
      !  Purpose: 
      !    This subroutine updates the density whenever called
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !


      subroutine updatedens(tfinal)
      use param
      implicit none
      
      integer :: ptype
      integer :: ix,iy,iz
      double precision, intent(in) :: tfinal
      
      ! retrieve particle type and subcell(s)
      ptype=nint(particle(7))
      ix=nint(particle(8))
      iy=nint(particle(9))
      iz=nint(particle(10))
      
      ! add the time the particle was in this subcell
      ncelldensnew3(ptype,ix,iy,iz)=
     &   ncelldensnew3(ptype,ix,iy,iz)+tfinal
      
      
      end subroutine updatedens
      

      
      