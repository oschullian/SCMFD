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


      subroutine readavp
      use param
      implicit none
      
      ! filename
      character(len=1024) :: filename
      ! format of the filename
      character(len=1024) :: format_string
      ! integer to test for reading 
      integer :: readerror
      
      ! loop integers
      integer :: ix, iy, iz
      integer :: id, ip
      
      ! read in integers
      integer :: ixr, iyr, izr, idr, ipr
      ! index in the array (1,noavp) and the read in one
      integer :: i, ir
      
      double precision :: avpveldummy 
      double precision, dimension(nvelchange) :: rvecd
      
      double precision :: x,y,chi,mu, meanv, ax,ay, an


      ! initialise value for avpvel if avpstart=0
      if (avpstart==0) then
        ! fill the whole cell with gas at T=tempwall
        avpvel3=0.d0
        ! loop over all subcells, dimensions and particle types
        do ix=1,ncellx
        do iy=1,ncelly
        do iz=1,ncellz
          do id=1,3
            do ip=1,nop
                
                ! call function that generates normal distributed random numbers in the 
                ! units of sigma
                !call normrnd(rvecd,nvelchange)

                ! multiply by sigma to get actual velocity and save
                avpvel3(1,ip,id,ix,iy,iz)=
     &        1.d0
              
            end do
          end do
        end do
        end do 
        end do 
        
        ! number of velocities that were saved in avpvel3
        nopcopos3=1
        
        ! set the arrays for updating the avpvel in the updating steps to 0
        avtpos3=0
        avtexchange3=0
        avpvelinter3=0.d0
        
      else 
        
        filename='avp.txt'
        ! open file
        open(unit=10, file=trim(filename), action='READ', status='old',
     & iostat=readerror)
      if (readerror.ne.0) then 
        write(*,*) 'can''t open avp file'
        write(*,*) 'check path'
        stop
      end if
      
      ! loop over all entries, check if they are correct and then save them
      do ip=1,nop
        do i=1,noavp
          do id=1,3
            do ix=1,ncellx
            do iy=1,ncelly
            do iz=1,ncellz
            
            read(10,*) ipr, ir, idr, ixr,iyr,izr, avpveldummy
            if (ipr.ne.ip) then
            write(*,*) 'read error'
            stop
            end if
            if (ir.ne.i) then
            write(*,*) 'read error'
            stop
            end if
            if (idr.ne.id) then
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
            avpvel3(i,ip,id,ix,iy,iz)=avpveldummy
            end do
            end do 
            end do
          end do
        end do
      end do
      
      ! close reading file
      close(unit=10)
     
      filename='ind.txt'
      ! open file with indexes (nopcopos) and read that in 
      open(unit=10, file=trim(filename), action='READ', status='old',
     & iostat=readerror)
      if (readerror.ne.0) then
        write(*,*) 'can''t ind file'
        write(*,*) 'check path'
        stop
      end if

      ! loop over all entries, check if they are correct and then save them
      do ip=1,nop
        do ix=1,ncellx
        do iy=1,ncelly
        do iz=1,ncellz
          read(10,*) ipr, ixr,iyr,izr, avpveldummy
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
          nopcopos3(ip,ix,iy,iz)=avpveldummy
        end do
        end do
        end do
      end do
      
      ! set the updating arrays to 0
      
      avtpos3=0
      avtexchange3=0
      avpvelinter3=0.d0
      ! close reading file
      close(unit=10)
      end if
      
      write(*,*) 'make done'

     
      end subroutine readavp
      

      
      
