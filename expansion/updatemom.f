      !
      !  Purpose: 
      !    This subroutine updates the vectors for the moments-calculation whenever called
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !


      subroutine updatemom(tfinal)
      use param
      implicit none
      double precision, intent(in) :: tfinal
      integer :: k,l
      integer :: count
      integer :: ptype
      integer :: ix,iy,iz
      
      ! retrieve type of the particle
      ptype=nint(particle(7))
      ix=nint(particle(8))
      iy=nint(particle(9))
      iz=nint(particle(10))
      
      ! the counter for the average (in future codes the density)
      n_count3(ptype,ix,iy,iz)=n_count3(ptype,ix,iy,iz)+tfinal
      
      ! first moment 
      do k=1,3
        mom_fi3(k,ptype,ix,iy,iz)=mom_fi3(k,ptype,ix,iy,iz)
     &   +particle(3+k)*tfinal
      end do 
      
      ! second moment
      count=0
      do k=1,3
        do l=k,3
          count=count+1
          mom_sec3(count,ptype,ix,iy,iz)=
     &      mom_sec3(count,ptype,ix,iy,iz)+
     &      particle(3+k)*particle(3+l)*tfinal;
        end do
      end do

      end subroutine updatemom
      

      
      
