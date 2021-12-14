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


      subroutine printavp
      use param
      implicit none
      
      character(len=1024) :: filename
      character(len=1024) :: format_string
      integer :: readerror
      !integer :: k,j,l,m, pt, nc, kf
      integer :: ip, ix, iy, iz, k, kf, id, i
      real, dimension(nvelchange) :: rvec
      integer :: indp, indexmax, indp2
      
      ! this is the counter for the filename, increase by one
      nprint=nprint+1
    
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
      
      ! reset the arrays for the exchange
      avtexchange3=0
      avtpos3=0
      
      
      
      ! this only generates the correct format for the filename, so that there 
      ! are no empty spaces in the filename
      if (nprint<10) then
        format_string="(A3,I1,A4)"
      else if (nprint<100) then
        format_string="(A3,I2,A4)"
      else if (nprint<1000) then
        format_string="(A3,I3,A4)"
      else if (nprint<10000) then
        format_string="(A3,I4,A4)"
      else if (nprint<100000) then
        format_string="(A3,I5,A4)"
      else if (nprint<1000000) then
        format_string="(A3,I6,A4)"
      end if
      
      ! generate the filename with the correct format
      write(filename,format_string) 'avp',nprint,'.txt'
      
      ! create or open the file for the values
      open(unit=10, file=trim(filename), action='write',
     & iostat=readerror)
      if (readerror.ne.0) then 
        write(*,*) 'can''t avp file to write'
        write(*,*) 'check path'
        stop
      end if
      
      
      ! save the values of the matrix in the format
      ! particletype, entry, direction(x,y,z), value
      do ip=1,nop
        do i=1,noavp
          do id=1,3
            do ix=1,ncellx
            do iy=1,ncelly
            do iz=1,ncellz
      write(10,*) ip, i, id, ix,iy,iz, avpvel3(i,ip,id,ix,iy,iz)
            end do 
            end do
            end do
          end do
        end do
      end do
      
      ! close the file
      close(unit=10)
      
      
      
      

      ! generate the filename with the correct format
      write(filename,format_string) 'ind',nprint,'.txt'

      ! create or open the file for the values
      open(unit=10, file=trim(filename), action='write',
     & iostat=readerror)
      if (readerror.ne.0) then
      write(*,*) 'can''t ind file to write'
      write(*,*) 'check path'
      stop
      end if


      ! save the values of the matrix in the format
      ! particletype, entry, direction(x,y,z), value
      do ip=1,nop
        do ix=1,ncellx
        do iy=1,ncelly
        do iz=1,ncellz
          write(10,*) ip, ix,iy,iz, nopcopos3(ip,ix,iy,iz)
        end do
      end do
      end do 
      end do
      ! close the file
      close(unit=10)
      write(*,*) 'print', nprint
      
      
     
     

      end subroutine printavp
      

      
      
