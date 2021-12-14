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
    
      !call newavp
      
      
      filename='avp.txt'
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
      
      
      
      filename='ind.txt'
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
      

      
      
