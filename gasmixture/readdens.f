      !
      !  Purpose: 
      !    This subroutine reads in the avp file in case one wants to propagate a simulation
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !


      subroutine readdens
      use param
      implicit none
      
      ! filename
      character(len=1024) :: filename
      ! format of the filename
      character(len=1024) :: format_string
      ! integer to test for reading 
      integer :: readerror
      
      ! loop integers (and the read in integers)
      integer :: ip, ix, iy, iz, ipr, ixr, iyr, izr
      ! read in density
      double precision :: densdummy 
      
      
      ! initialise value for avpvel if avpstart=0
      if (avpstart==0) then
        ! in this case we fill the whole cell with equal density 
        ! set all values to zero and start
        do ip=1,nop
          ncelldens3(ip,:,:,:)=npartdens(ip)
        end do
        ! initialise the update dens to 0
        ncelldensnew3=0.d0
        
      else 
        ! define format of filename (so no empty spaces are there)
        if (avpstart<10) then
          format_string="(A4,I1,A4)"
        else if (avpstart<100) then
          format_string="(A4,I2,A4)"
        else if (avpstart<1000) then
          format_string="(A4,I3,A4)"
        else if (avpstart<10000) then
          format_string="(A4,I4,A4)"
        else if (avpstart<100000) then
          format_string="(A4,I5,A4)"
        else if (avpstart<1000000) then
          format_string="(A4,I6,A4)"
        end if
        
        ! create filename
        write(filename,format_string) 'dens',avpstart,'.txt'
        
        ! open file
        open(unit=10, file=trim(filename), action='READ', status='old',
     & iostat=readerror)
      if (readerror.ne.0) then 
        write(*,*) 'can''t open dens file'
        write(*,*) 'check path'
        stop
      end if
      
      ! loop over all entries, check if they are correct and then save them
      do ip=1,nop
            do ix=1,ncellx
            do iy=1,ncelly
            do iz=1,ncellz
            read(10,*) ipr, ixr, iyr, izr, densdummy
            if (ipr.ne.ip) then
            write(*,*) 'read error'
            stop
            end if
            if (ixr.ne.ix) then
            write(*,*) 'read error'
            stop
            end if
            if (iyr.ne.iy) then
            write(*,*) 'read error'
            stop
            end if
            if (izr.ne.iz) then
            write(*,*) 'read error'
            stop
            end if
            ncelldens3(ip,ix,iy,iz)=densdummy
            end do
            end do 
            end do
      end do
      
      ncelldensnew3=0.d0
      
      
      ! close reading file
      close(unit=10)

      end if

      end subroutine readdens
      

      
      