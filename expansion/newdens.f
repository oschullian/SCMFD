      !
      !  Purpose: 
      !    This subroutine prints out all the the avp, for investigation
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !


      subroutine newdens
      use param
      implicit none
      
      character(len=1024) :: filename
      character(len=1024) :: format_string
      integer :: readerror
      integer :: ip, ix, iy, iz
      double precision :: volsubcell
      


      volsubcell=product(distall)/dble(product(ncellall))
      

      do ip=1,nop
        ncelldensnew3(ip,:,:,:)=ncelldensnew3(ip,:,:,:)*
     &   influx(ip)/(volsubcell*dble(np(ip)*modnewavp))
      end do

      
      !adjust the density
      do ip=1,nop
      ncelldens3(ip,:,:,:)=ncelldens3(ip,:,:,:)*(1.d0-inpcellperc)+
     &     ncelldensnew3(ip,:,:,:)*inpcellperc
      
      end do
      ! set the density for update to 0
      ncelldensnew3=0.d0
      
      
      end subroutine newdens
      

      
      
