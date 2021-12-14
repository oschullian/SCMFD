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


      subroutine newavp
      use param
      implicit none
      
      character(len=1024) :: filename
      character(len=1024) :: format_string
      integer :: readerror
      !integer :: k,j,l,m, pt, nc, kf
      integer :: ip, ix, iy, iz, k, kf, id, i
      real, dimension(nvelchange) :: rvec
      integer :: indp, indexmax, indp2
      double precision, dimension(nop) :: propexchange

      ! move the velocities in one column by avtexchange position downwards
      ! discard the ones that would lie outside the array
      do ip=1,nop
      do ix=1,ncellx
      do iy=1,ncelly
      do iz=1,ncellz
      do k=noavp-avtexchange3(ip,ix,iy,iz),1,-1
      kf=k+avtexchange3(ip,ix,iy,iz)
      avpvel3(kf,ip,:,ix,iy,iz)=avpvel3(k,ip,:,ix,iy,iz)
      end do
      ! save the avtexchange velocities in the avpvelinter onto the new "free" 
      ! positions in avpvel3
      do k=avtexchange3(ip,ix,iy,iz),1,-1
      kf=avtexchange3(ip,ix,iy,iz)-k+1
      avpvel3(kf,ip,:,ix,iy,iz)=avpvelinter3(k,ip,:,ix,iy,iz)
      end do
      
      ! adjust the nopcopos array 
      nopcopos3(ip,ix,iy,iz)=
     &  min(nopcopos3(ip,ix,iy,iz)+avtexchange3(ip,ix,iy,iz),noavp)
      end do
      end do
      end do 
      end do
      
      do ip=1,nop
      propexchange(ip)=sum(avtexchange3(ip,:,:,:))/
     &(product(ncellall)*inpcellperc*noavp)
      write(122,*) propexchange
      end do
 
      ! reset the arrays for the exchange
      avtexchange3=0
      avtpos3=0
      
      
      
      
      
     
     

      end subroutine newavp
      

      
      
