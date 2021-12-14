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


      subroutine printdens
      use param
      implicit none
      
      character(len=1024) :: filename
      character(len=1024) :: format_string
      integer :: readerror
      integer :: ip, ix, iy, iz
      
      ! this is the counter for the filename, increase by one
      nprintdens=nprintdens+1
      
      ! generate the filename with the correct format
      filename='dens.txt'
      ! create or open the file for the values
      open(unit=10, file=trim(filename), action='write',
     & iostat=readerror)
      if (readerror.ne.0) then 
        write(*,*) 'can''t dens file to write'
        write(*,*) 'check path'
        stop
      end if
     
      ! adjust the density 
      do ip=1,nop
        ncelldens3(ip,:,:,:)=ncelldens3(ip,:,:,:)*(1.d0-inpcellperc)+
     &     ncelldensnew3(ip,:,:,:)*
     &     dble(ncellx*ncelly*ncellz)/sum(ncelldensnew3(ip,:,:,:))*
     &     npartdens(ip)*inpcellperc
     
      end do 
      ! set the density for update to 0
      ncelldensnew3=0.d0
      
      
      ! save the values of the matrix in the format
      ! particletype, entry, direction(x,y,z), value
      do ip=1,nop
            do ix=1,ncellx
            do iy=1,ncelly
            do iz=1,ncellz
              write(10,*) ip, ix, iy, iz, 
     &    ncelldens3(ip,ix, iy, iz)
            end do
            end do
            end do
      end do
      
      ! close the file
      close(unit=10)

      ! print out the number of nprint so one can check the newest (altough ls -lrt will
      ! also do
      write(*,*) 'printdens', nprintdens
      
      end subroutine printdens
      

      
      
